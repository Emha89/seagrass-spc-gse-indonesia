/************************************************************
 * 03_GEE_app_viewer.js
 *
 * Purpose:
 *   Interactive web application for visual exploration of
 *   two-stage GSE-based seagrass percent cover (SPC)
 *   predictions across five study sites in Indonesia.
 *
 *   Features:
 *     - Select study site and survey year
 *     - Toggle training point overlay
 *     - Click map to query predicted SPC value
 *     - SPC colour legend
 *
 *   The app is for qualitative visual assessment only,
 *   as described in Section 3.5 of the manuscript.
 *
 *   Live app:
 *   https://muhammadhafizt.users.earthengine.app/view/
 *   seagrasspercentcoverindonesia
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

// SPC prediction assets (output from 02_SPC_spatial_deployment.js)
var SPC_ASSETS = {
  '2018': 'projects/YOUR_PROJECT/assets/SPC_GSE_twoStage_2018',
  '2019': 'projects/YOUR_PROJECT/assets/SPC_GSE_twoStage_2019',
  '2021': 'projects/YOUR_PROJECT/assets/SPC_GSE_twoStage_2021',
  '2022': 'projects/YOUR_PROJECT/assets/SPC_GSE_twoStage_2022'
};

// Study area boundaries (FeatureCollection with 'loc' property)
var AOI_ASSET = 'projects/YOUR_PROJECT/assets/YOUR_STUDY_AREA';

// Field survey training points
var TRAINING_PTS_ASSET = 'projects/YOUR_PROJECT/assets/YOUR_TRAINING_POINTS';

// ==========================================================
// LOAD ASSETS
// ==========================================================

var AOI_FC       = ee.FeatureCollection(AOI_ASSET);
var TRAINING_PTS = ee.FeatureCollection(TRAINING_PTS_ASSET);

// ==========================================================
// SITE CENTRES [longitude, latitude, zoom]
// ==========================================================

var SITES = {
  'Full Indonesia': [118.0,   -2.5,   5],
  'Ayau':           [131.047,  0.362, 15],
  'Bintan':         [104.569,  1.232, 15],
  'Karimunjawa':    [110.479, -5.772, 15],
  'Komodo':         [119.726, -8.568, 15],
  'Rote':           [122.811,-10.778, 15]
};

var YEARS = ['2018', '2019', '2021', '2022'];

// ==========================================================
// VISUALISATION PARAMETERS
// ==========================================================

var SPC_VIS = {
  min:     0,
  max:     100,
  palette: ['ffffff', 'e5f5e0', 'a1d99b', '31a354', '006d2c']
};

// ==========================================================
// MAP WIDGET
// ==========================================================

var mapWidget = ui.Map();
mapWidget.setOptions('SATELLITE');
mapWidget.setCenter(118, -2, 5);
mapWidget.style().set('cursor', 'crosshair');

var clickLabel = ui.Label({
  value: '',
  style: {
    fontSize:        '13px',
    fontWeight:      'bold',
    backgroundColor: 'white',
    padding:         '4px',
    margin:          '4px'
  }
});
mapWidget.add(clickLabel);

// ==========================================================
// UI ELEMENTS
// ==========================================================

var title = ui.Label({
  value: 'Seagrass Percent Cover (Two-stage RF, GSE-based)',
  style: { fontSize: '15px', fontWeight: 'bold', margin: '0 0 4px 0' }
});

var subtitle = ui.Label({
  value: 'Interactive visualisation for qualitative assessment only',
  style: { fontSize: '11px', color: '#555555', margin: '0 0 10px 0' }
});

var locationSelect = ui.Select({
  items:       Object.keys(SITES),
  placeholder: 'Select location'
});

var yearSelect = ui.Select({
  items: YEARS,
  value: YEARS[0]
});

var trainingToggle = ui.Checkbox({
  label: 'Show training points',
  value: false
});

var infoLabel = ui.Label({
  value: 'Select a location to display SPC.',
  style: { fontSize: '12px', color: '#aa0000', margin: '6px 0 0 0' }
});

// ==========================================================
// LEGEND
// ==========================================================

var legend = ui.Panel({ style: { margin: '12px 0 0 0' } });
legend.add(ui.Label({ value: 'SPC (%)', style: { fontWeight: 'bold', margin: '0 0 4px 0' } }));

var palette = SPC_VIS.palette;
var labels  = ['0', '25', '50', '75', '100'];

for (var i = 0; i < palette.length; i++) {
  legend.add(ui.Panel({
    widgets: [
      ui.Label('', {
        backgroundColor: palette[i],
        padding: '8px',
        margin:  '2px 6px 2px 0'
      }),
      ui.Label(labels[i], { margin: '4px 0' })
    ],
    layout: ui.Panel.Layout.Flow('horizontal')
  }));
}

// ==========================================================
// MAP UPDATE FUNCTION
// ==========================================================

var lastCenteredLoc = null;

function updateMap() {
  mapWidget.layers().reset();

  var loc  = locationSelect.getValue();
  var year = yearSelect.getValue();

  if (!loc) {
    infoLabel.setValue('Select a location to display SPC.');
    return;
  }

  infoLabel.setValue('Loading...');

  if (loc !== lastCenteredLoc) {
    var siteInfo = SITES[loc];
    if (siteInfo) mapWidget.setCenter(siteInfo[0], siteInfo[1], siteInfo[2]);
    lastCenteredLoc = loc;
  }

  if (loc === 'Full Indonesia') {
    mapWidget.addLayer(
      AOI_FC.style({ color: '000000', fillColor: '00000000', width: 2 }),
      {}, 'AOI'
    );
    infoLabel.setValue('');
    return;
  }

  var aoi = AOI_FC.filter(ee.Filter.eq('loc', loc));
  var img = ee.Image(SPC_ASSETS[year]).clip(aoi);

  mapWidget.addLayer(img, SPC_VIS, 'SPC ' + year);
  mapWidget.addLayer(
    aoi.style({ color: '000000', fillColor: '00000000', width: 2 }),
    {}, 'AOI'
  );

  if (trainingToggle.getValue()) {
    var pts = TRAINING_PTS.filter(ee.Filter.eq('loc', loc));
    mapWidget.addLayer(
      pts.style({ color: 'ff0000', pointSize: 3 }),
      {}, 'Training points'
    );
  }

  infoLabel.setValue('');
}

// ==========================================================
// CLICK INSPECTOR
// ==========================================================

mapWidget.onClick(function(coords) {
  var loc  = locationSelect.getValue();
  var year = yearSelect.getValue();
  if (!loc || loc === 'Full Indonesia') return;

  clickLabel.setValue('SPC: loading...');

  var img   = ee.Image(SPC_ASSETS[year]);
  var point = ee.Geometry.Point([coords.lon, coords.lat]);

  var value = img.reduceRegion({
    reducer:    ee.Reducer.first(),
    geometry:   point,
    scale:      10,
    bestEffort: true
  });

  value.evaluate(function(v) {
    if (v && v.SPC_pred !== null && v.SPC_pred !== undefined) {
      clickLabel.setValue('SPC = ' + v.SPC_pred.toFixed(1) + ' %');
    } else {
      clickLabel.setValue('No data at this point');
    }
  });

  mapWidget.layers().forEach(function(l) {
    if (l.getName() === 'Click point') mapWidget.layers().remove(l);
  });
  mapWidget.addLayer(point, { color: 'FF0000' }, 'Click point');
});

// ==========================================================
// WIRE UI EVENTS
// ==========================================================

locationSelect.onChange(updateMap);
yearSelect.onChange(updateMap);
trainingToggle.onChange(updateMap);

// ==========================================================
// LAYOUT
// ==========================================================

var panel = ui.Panel({
  widgets: [
    title,
    subtitle,
    ui.Label('Location', { fontWeight: 'bold' }),
    locationSelect,
    ui.Label('Year', { fontWeight: 'bold', margin: '8px 0 0 0' }),
    yearSelect,
    ui.Label('', { margin: '4px 0' }),
    trainingToggle,
    infoLabel,
    legend
  ],
  style: { width: '300px', padding: '12px' }
});

ui.root.clear();
ui.root.add(ui.SplitPanel({
  firstPanel:  panel,
  secondPanel: mapWidget,
  orientation: 'horizontal',
  wipe:        false
}));
