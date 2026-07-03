/************************************************************
 * 02_SPC_spatial_deployment.js
 *
 * Purpose:
 *   Spatial deployment of the two-stage Random Forest SPC
 *   prediction model using Google Satellite Embeddings (GSE).
 *
 *   Stage A: RF morphology classification (3 classes) using
 *            GSE dimensions. Class membership probabilities
 *            are extracted per pixel.
 *   Stage B: RF SPC regression using GSE dimensions combined
 *            with Stage A morphology class probabilities.
 *
 *   Predictions are masked to the seagrass persistence mask
 *   and exported to GEE Assets and Google Drive for each
 *   survey year (2018, 2019, 2021, 2022).
 *
 *   Outputs correspond to Figure 8 in the manuscript.
 *   This script is for spatial plausibility assessment only
 *   and is aligned conceptually with the R-based LOLO
 *   cross-validation (STEP_05_two_stage_RF_LOLOCV.R).
 *
 * Study:
 *   Interpretable Estimation of Seagrass Percent Cover Across
 *   Extensive Coastal Environments in Indonesia Using Google
 *   Satellite Embeddings
 *   International Journal of Digital Earth (2026)
 *   Manuscript ID: TJDE-2026-0306
 *
 * Authors: Muhammad Hafizt, Stuart Phinn, Pramaditya Wicaksono,
 *          Udhi Eko Hernawan, Huwaida Nur Salsabila,
 *          Mitchell Lyons, Kathryn McMahon, Chris Roelfsema
 *
 * NOTE:
 *   This script requires access to GEE assets that are not
 *   publicly available. See README_GEE.md for the full list
 *   of required assets and how to obtain or substitute them.
 *   Update the USER CONFIGURATION section before running.
 ************************************************************/

// ==========================================================
// USER CONFIGURATION
// ==========================================================

// Field survey training points with morph3 and total_SPC properties
// Data available from corresponding author upon request (see paper)
var TRAIN_PTS_ASSET = 'projects/YOUR_PROJECT/assets/YOUR_TRAINING_POINTS';

// Study area boundary (FeatureCollection)
var AOI_ASSET = 'projects/YOUR_PROJECT/assets/YOUR_STUDY_AREA';

// Seagrass persistence mask (binary raster: 1 = seagrass present)
// Generated from multi-year seagrass occurrence mapping
var PERSISTENCE_MASK_ASSET = 'projects/YOUR_PROJECT/assets/YOUR_PERSISTENCE_MASK';

// GEE Asset path prefix for exported SPC prediction rasters
var EXPORT_ASSET_PREFIX = 'projects/YOUR_PROJECT/assets/YOUR_OUTPUT_PREFIX_';

// Google Drive export folder
var DRIVE_FOLDER = 'your_export_folder';

// Survey years to process
var YEARS = [2018, 2019, 2021, 2022];

// Pixel scale (metres)
var SCALE = 10;

// Random seed (matched to R scripts)
var SEED = 42;

// ==========================================================
// RF HYPERPARAMETERS
// (matched conceptually to R grid search results)
// ==========================================================

var RF_CLS = {
  numberOfTrees:      700,
  variablesPerSplit:  10,
  minLeafPopulation:  3,
  bagFraction:        0.7,
  seed:               SEED
};

var RF_REG = {
  numberOfTrees:      700,
  variablesPerSplit:  10,
  minLeafPopulation:  5,
  bagFraction:        0.7,
  seed:               SEED
};

// ==========================================================
// LOAD ASSETS
// ==========================================================

var AOI             = ee.FeatureCollection(AOI_ASSET).geometry();
var PERSISTENCE_MASK = ee.Image(PERSISTENCE_MASK_ASSET);

// GSE band names (A00-A63)
var GSE_BANDS = ee.List.sequence(0, 63).map(function(i) {
  return ee.String('A').cat(ee.Number(i).format('%02d'));
});

// Morphology class codes (must match R script reclassification):
//   0 = mixed_short_plus_mono_short
//   1 = mixed_long
//   2 = mono_Ea
var MORPH_CLASSES = [
  'mixed_short_plus_mono_short',
  'mixed_long',
  'mono_Ea'
];

// ==========================================================
// LOAD AND PREPARE TRAINING DATA
// ==========================================================

var rawPts = ee.FeatureCollection(TRAIN_PTS_ASSET)
  .filter(ee.Filter.gte('total_SPC', 0))
  .filter(ee.Filter.lte('total_SPC', 100));

// Reclassify morphology to numeric codes (matching R morph3 reclassification)
function addMorph3Code(f) {
  var sg = ee.String(f.get('sg_morpho'));
  var ea = ee.Number(f.get('Ea_SPC'));

  var morph3 = ee.Algorithms.If(
    sg.equals('mixed_long'), 'mixed_long',
    ee.Algorithms.If(
      sg.equals('mono'),
      ee.Algorithms.If(ea.gt(0), 'mono_Ea', 'mixed_short_plus_mono_short'),
      ee.Algorithms.If(sg.equals('mixed_short'), 'mixed_short_plus_mono_short', null)
    )
  );

  var morph3Code = ee.Algorithms.If(
    ee.String(morph3).equals('mixed_short_plus_mono_short'), 0,
    ee.Algorithms.If(
      ee.String(morph3).equals('mixed_long'), 1,
      ee.Algorithms.If(ee.String(morph3).equals('mono_Ea'), 2, null)
    )
  );

  return f.set({ morph3: morph3, morph3_code: morph3Code });
}

// ==========================================================
// GSE LOADING FUNCTION
// ==========================================================

function loadGSE(year) {
  return ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
    .filterDate(
      ee.Date.fromYMD(year, 1, 1),
      ee.Date.fromYMD(year, 1, 1).advance(1, 'year'))
    .filterBounds(AOI)
    .mosaic()
    .select(GSE_BANDS);
}

// ==========================================================
// SAMPLE GSE AT TRAINING POINTS PER YEAR
// ==========================================================

function sampleYear(year) {
  var img = loadGSE(year);
  var pts = rawPts.filter(ee.Filter.eq('year', year));
  return img.sampleRegions({
    collection: pts,
    properties: pts.first().propertyNames(),
    scale:      SCALE,
    geometries: true
  });
}

// Build multi-year training table with morph3 reclassification
var training = ee.FeatureCollection(YEARS.map(sampleYear)).flatten()
  .map(addMorph3Code)
  .filter(ee.Filter.notNull(['morph3_code']));

print('Training samples after GSE extraction and morph3 reclass:', training.size());
print('Morph3 distribution:', training.aggregate_histogram('morph3'));

// ==========================================================
// STAGE A: MORPHOLOGY CLASSIFICATION RF
// ==========================================================

var rfMorph = ee.Classifier.smileRandomForest(RF_CLS)
  .train({
    features:        training,
    classProperty:   'morph3_code',
    inputProperties: GSE_BANDS
  });

print('Morphology RF trained');

// Add per-class probabilities to training set (MULTIPROBABILITY mode)
var trainingWithProbArr = training.classify(
  rfMorph.setOutputMode('MULTIPROBABILITY'), 'prob_arr'
);

var trainingWithProb = trainingWithProbArr.map(function(f) {
  var arr = ee.Array(f.get('prob_arr'));
  return f.set({
    prob_mixed_short_plus_mono_short: arr.get([0]),
    prob_mixed_long:                  arr.get([1]),
    prob_mono_Ea:                     arr.get([2])
  });
});

// ==========================================================
// STAGE B: SPC REGRESSION RF (GSE + morphology probabilities)
// ==========================================================

var REG_BANDS = GSE_BANDS.cat([
  'prob_mixed_short_plus_mono_short',
  'prob_mixed_long',
  'prob_mono_Ea'
]);

var rfSPC = ee.Classifier.smileRandomForest(RF_REG)
  .setOutputMode('REGRESSION')
  .train({
    features:        trainingWithProb,
    classProperty:   'total_SPC',
    inputProperties: REG_BANDS
  });

print('SPC regression RF trained');

// ==========================================================
// APPLY MODEL PER YEAR AND EXPORT
// ==========================================================

YEARS.forEach(function(y) {

  print('Applying model for year:', y);

  var gse = loadGSE(y);

  // Stage A: morphology class probabilities
  var probImg = gse.classify(
    rfMorph.setOutputMode('MULTIPROBABILITY')
  ).arrayFlatten([[
    'prob_mixed_short_plus_mono_short',
    'prob_mixed_long',
    'prob_mono_Ea'
  ]]);

  // Stage B: SPC regression
  var spcPred = gse
    .addBands(probImg)
    .select(REG_BANDS)
    .classify(rfSPC)
    .rename('SPC_pred');

  // Apply persistence mask and clip to AOI (for display and export)
  var spcMasked = spcPred
    .updateMask(PERSISTENCE_MASK)
    .clip(AOI);

  // Map display
  Map.addLayer(
    spcMasked,
    {
      min:     0,
      max:     100,
      palette: ['ffffff', 'e5f5e0', 'a1d99b', '31a354', '006d2c']
    },
    'Predicted SPC ' + y,
    false
  );

  // Export to GEE Asset
  Export.image.toAsset({
    image:       spcMasked,
    description: 'SPC_GSE_twoStage_' + y,
    assetId:     EXPORT_ASSET_PREFIX + y,
    region:      AOI,
    scale:       SCALE,
    maxPixels:   1e13
  });

  // Export to Google Drive (for QGIS / local use)
  Export.image.toDrive({
    image:           spcMasked,
    description:     'SPC_GSE_twoStage_' + y + '_Drive',
    folder:          DRIVE_FOLDER,
    fileNamePrefix:  'SPC_GSE_twoStage_' + y,
    region:          AOI,
    scale:           SCALE,
    maxPixels:       1e13
  });

  print('Export queued (Asset + Drive) for year:', y);
});

Map.centerObject(AOI, 6);
