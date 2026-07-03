/************************************************************
 * 01_S2_GSE_extraction.js
 *
 * Purpose:
 *   Extract Sentinel-2 percentile composites and Google Satellite
 *   Embedding (GSE) values at field survey locations, then export
 *   as CSV for use in R-based sensitivity analysis and modelling
 *   (STEP_01 through STEP_05 R scripts).
 *
 *   Three S2 stacks are generated per year:
 *     - raw:  no sunglint correction
 *     - post: sunglint correction applied after stacking
 *     - pre:  sunglint correction applied before stacking
 *
 *   GSE is extracted for the year matching the field survey year
 *   and joined to S2 samples by unique observation ID (gee_id).
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

// Field survey training points (FeatureCollection with gee_id property)
// Data available from corresponding author upon request (see paper)
var TRAIN_PTS_ASSET = 'projects/YOUR_PROJECT/assets/YOUR_TRAINING_POINTS';

// Study area boundary (FeatureCollection or Geometry)
var TRAINING_EXTENT_ASSET = 'projects/YOUR_PROJECT/assets/YOUR_STUDY_AREA';

// Bathymetry raster (used for depth stratification in spectral profiles)
var DEPTH_ASSET = 'projects/YOUR_PROJECT/assets/YOUR_BATHYMETRY';

// Distance to land raster
var DIST_TO_LAND_ASSET = 'projects/coral_atlas/global_datasets/osm_distToLand_indo';

// GSE year to extract (match to field survey year)
var GSE_YEAR = 2019;

// Sentinel-2 composite date range (match to field survey year)
var S2_START = '2019-01-01';
var S2_END   = '2019-12-31';

// Cloud cover threshold for S2 scene filtering (%)
var CLOUDY_THRESHOLD = 60;

// Google Drive export folder
var DRIVE_FOLDER = 'your_export_folder';

// Pixel scale (metres)
var SCALE = 10;

// ==========================================================
// LOAD ASSETS
// ==========================================================

var trainingPts     = ee.FeatureCollection(TRAIN_PTS_ASSET);
var trainingExtent  = ee.FeatureCollection(TRAINING_EXTENT_ASSET).geometry();
var depth           = ee.Image(DEPTH_ASSET).multiply(-1).rename('depth');
var distToLand      = ee.Image(DIST_TO_LAND_ASSET).rename('distToLand');

print('Training points:', trainingPts.size());

// ==========================================================
// S2 BANDS TO KEEP
// ==========================================================

var S2_KEEP_BANDS = ['B2', 'B3', 'B4', 'B8'];

// Training extent mask (used to restrict cloud masking to study area)
var trainingExtentMask = ee.Image().byte()
  .paint(ee.Feature(trainingExtent, {zone: 1}), 'zone');

// ==========================================================
// CLOUD MASKING FUNCTION (Sentinel-2 Level-2A SCL)
// ==========================================================

function maskS2clouds(keepBands, extentMask) {
  return function(image) {
    var scl        = image.select('SCL');
    var cloudMask  = scl.lt(7).and(scl.gt(3));
    return image.select(keepBands)
      .updateMask(cloudMask)
      .updateMask(extentMask);
  };
}

// ==========================================================
// SUNGLINT CORRECTION FUNCTIONS
// ==========================================================

// Applied after stacking (to percentile composite)
function applySunglintPost(image) {
  var corrected = image.select(['B4_p60', 'B3_p60', 'B2_p60'])
    .subtract(image.select('B8_p60'))
    .rename(['B4_p60', 'B3_p60', 'B2_p60']);
  return image.addBands(corrected, null, true);
}

// Applied before stacking (to each individual S2 image)
function applySunglintPre(image) {
  var corrected = image.select(['B4', 'B3', 'B2'])
    .subtract(image.select('B8'))
    .rename(['B4', 'B3', 'B2']);
  return image.addBands(corrected, null, true);
}

// ==========================================================
// BUILD SENTINEL-2 STACKS
// ==========================================================

var s2Collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(trainingExtent)
  .filterDate(S2_START, S2_END)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUDY_THRESHOLD))
  .map(maskS2clouds(S2_KEEP_BANDS, trainingExtentMask));

// Standard deviation bands (temporal variability within year)
var s2StdDev = s2Collection.reduce(ee.Reducer.stdDev())
  .rename(S2_KEEP_BANDS.map(function(b) {
    return ee.String(b).cat('_stddev');
  }));

// [1] Raw stack — no sunglint correction
var s2Raw = s2Collection.reduce(ee.Reducer.median()).uint16()
  .addBands(s2Collection.reduce(
    ee.Reducer.percentile(ee.List.sequence(0, 100, 20))).uint16())
  .addBands(s2StdDev);

// [2] Post-stack sunglint correction
var s2Post = applySunglintPost(s2Raw);

// [3] Pre-stack sunglint correction
var s2CollectionCorrected = s2Collection.map(applySunglintPre);
var s2Pre = s2CollectionCorrected.reduce(ee.Reducer.median()).uint16()
  .addBands(s2CollectionCorrected.reduce(
    ee.Reducer.percentile(ee.List.sequence(0, 100, 20))).uint16())
  .addBands(s2StdDev);

// ==========================================================
// ADD DERIVED FEATURES
// ==========================================================

function addDerivedFeatures(stack) {
  var slope    = ee.Terrain.slope(depth).rename('slope');
  var rugosity = slope.reduceNeighborhood({
    reducer: ee.Reducer.stdDev(),
    kernel:  ee.Kernel.circle(100, 'meters')
  }).unmask(0).rename('rugosity');

  return stack
    .addBands(stack.select('B4_p60').divide(stack.select('B3_p60')).rename('rg_median'))
    .addBands(stack.select('B4_p60').divide(stack.select('B2_p60')).rename('rb_median'))
    .addBands(stack.normalizedDifference(['B8_p60', 'B4_p60']).rename('ndvi'))
    .addBands(stack.normalizedDifference(['B3_p60', 'B8_p60']).rename('ndwi1'))
    .addBands(stack.normalizedDifference(['B4_p60', 'B3_p60']).rename('ndti'))
    .addBands(depth)
    .addBands(slope)
    .addBands(rugosity)
    .addBands(distToLand);
}

var s2RawEnh  = addDerivedFeatures(s2Raw);
var s2PostEnh = addDerivedFeatures(s2Post);
var s2PreEnh  = addDerivedFeatures(s2Pre);

// ==========================================================
// SAMPLE S2 AT TRAINING POINTS
// ==========================================================

function sampleS2(stack, label) {
  var KEY = 'gee_id';
  var sampled = stack.sampleRegions({
    collection: trainingPts,
    properties: [KEY],
    scale:      SCALE,
    geometries: true
  });

  var withCoords = sampled.map(function(f) {
    var coords = f.geometry().centroid(10).coordinates();
    return f.set({
      lon: ee.List(coords).get(0),
      lat: ee.List(coords).get(1)
    });
  });

  var clean = withCoords
    .filter(ee.Filter.notNull(['lon', 'lat']))
    .map(function(f) { return f.set(KEY, ee.String(f.get(KEY))); });

  print('S2 samples (' + label + '):', clean.size());
  return clean;
}

var s2SamplesRaw  = sampleS2(s2RawEnh,  's2_raw');
var s2SamplesPost = sampleS2(s2PostEnh, 's2_post');
var s2SamplesPre  = sampleS2(s2PreEnh,  's2_pre');

// ==========================================================
// EXTRACT GSE AT TRAINING POINTS
// ==========================================================

function buildGSE(year) {
  var img = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
    .filterDate(
      ee.Date.fromYMD(year, 1, 1),
      ee.Date.fromYMD(year + 1, 1, 1))
    .mosaic()
    .unmask(-9999, false);
  return img.rename(img.bandNames().map(function(b) {
    return ee.String('GSE_').cat(b);
  }));
}

var gseImg = buildGSE(GSE_YEAR)
  .reproject({ crs: s2Raw.projection(), scale: SCALE });

print('GSE year:', GSE_YEAR);

var gseSamples = gseImg.sampleRegions({
  collection: trainingPts.map(function(f) {
    return f.set('gee_id', ee.String(f.get('gee_id')));
  }),
  properties: ['gee_id'],
  scale:      SCALE,
  geometries: false
}).map(function(f) {
  return f.set('gee_id', ee.String(f.get('gee_id')));
});

print('GSE samples:', gseSamples.size());

// ==========================================================
// JOIN S2 AND GSE BY gee_id
// ==========================================================

function joinS2andGSE(s2Samples, gseImage, label) {
  var KEY     = 'gee_id';
  var gseCols = gseImage.bandNames();

  var saveFirst = ee.Join.saveFirst('match');
  var filter    = ee.Filter.equals({ leftField: KEY, rightField: KEY });
  var joined    = ee.FeatureCollection(saveFirst.apply(s2Samples, gseSamples, filter));

  var fallbackVals = gseCols.map(function(_) { return -9999; });
  var fallback     = ee.Feature(null, ee.Dictionary.fromLists(gseCols, fallbackVals));

  var merged = joined.map(function(f) {
    var match = ee.Feature(f.get('match'));
    var src   = ee.Algorithms.If(match, match, fallback);
    return f.copyProperties(ee.Feature(src), gseCols);
  });

  print('[' + label + '] merged size:', merged.size());

  // Export as CSV (primary input for R scripts)
  Export.table.toDrive({
    collection:     merged,
    description:    label + '_S2_GSE_joined',
    folder:         DRIVE_FOLDER,
    fileNamePrefix: label + '_S2_GSE_joined',
    fileFormat:     'CSV'
  });
}

joinS2andGSE(s2SamplesRaw,  gseImg, 's2_raw');
joinS2andGSE(s2SamplesPost, gseImg, 's2_post');
joinS2andGSE(s2SamplesPre,  gseImg, 's2_pre');

// ==========================================================
// MAP DISPLAY
// ==========================================================

Map.centerObject(trainingExtent, 6);
Map.addLayer(s2SamplesRaw, { color: 'red' },   'Samples - raw');
Map.addLayer(s2SamplesPost, { color: 'green' }, 'Samples - post');
Map.addLayer(s2SamplesPre, { color: 'blue' },  'Samples - pre');
