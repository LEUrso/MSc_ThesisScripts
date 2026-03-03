/**** ================================ LAST WORKFLOW UPDATE: JAN 19, 2026
 * Landsat image preprocesser version 2.2
 * Used to preprocess validation images 
 * incorporates updates to training image preprocessor
 *  
 * loads Landsat C2 T1_TOA scenes from SCENE_IDS
 * applies cloud mask 
 * applies rock mask (from Moussavi et al. 2022)
 * clips to AOI
 * adds NDWI_ICE band = normDiff(B2,B4) (improvment from v1)
 * exports preprocessed scenes (Assets or Drive) * 

/**** version info ===========================
 * training image preprocesser version 2
 * create:           14 jan 2026
 * ---
 * this script was created by Luke Urso as part of MSc thesis 
 ===============================================****/
 
// 1) SCENE LIST  (one blockat a time)
var SCENE_IDS = [

//OST HARDER validation scenes
  'LC08_033001_20210730',
  'LC08_033246_20210730',
  'LC08_033001_20180807',
  'LC08_035246_20180805',
  'LC08_039245_20180801',
  'LC08_023248_20160710',
  'LC08_047244_20220719',
  'LC08_033246_20220717',

//PTM secondary validation scenes
//  'LC09_036002_20250807',
//  'LC09_066242_20240619',
//  'LC08_036002_20190714',
//  'LC08_036002_20170622',
//  'LC08_036002_20150617',
//  'LC08_037002_20210726',

];


// 2) AOI asset 
// OSTENFELD
var AOI_ASSET = 'projects/vernal-signal-270100/assets/StudyArea/Mannual_RockMask/OST_STUDYAREA_ROCKMASK_SIMPLE';
// PETERMANN
//var AOI_ASSET = 'projects/vernal-signal-270100/assets/StudyArea/Mannual_RockMask/PTM_LONGAOI_WFJORDC_ROCKMASK';

// 3) Export destination + mode
// 'DRIVE' or 'ASSET'
var EXPORT_MODE = 'ASSET';

// Drive settings (only used if EXPORT_MODE === 'DRIVE')
var DRIVE_FOLDER = 'TrainingImages';

// Asset settings (only used if EXPORT_MODE === 'ASSET')
// create the destination folder/collection path in Assets first.
var DEST_ASSET_PREFIX =
  'projects/vernal-signal-270100/assets/SceneAssets/ValidationScenes/HarderErrorScenes_OST/';

// Name prefix for exported images/tasks (edit to OST_/PTM_ as needed)
var NAME_PREFIX = 'OST_PVI_'; // e.g., 'OST_PVI_' or 'PTM_PVI_'

// 4) Export parameters
// ------------------------
var EXPORT_SCALE_M = 30;
var MAX_PIXELS     = 1e13;

// 5) Optional visualization
// ------------------------
var ADD_RGB_LAYERS = true;

var VIZ_RGB = { bands: ['B4', 'B3', 'B2'], min: 0, max: 1.0 };

// NDWI_ICE viz (adjust as needed)
var VIZ_NDWI_ICE = {
  min: 0.115,
  max: 0.22,
  palette: ['181a1b','697578','00ffe7']
};

/*******************************************************
 * MODULES
 *******************************************************/
var cloud_mask = require('users/LukeUrso/GEEScripts:ProcessingModules/LC08_CloudMask');
var rock_mask  = require('users/LukeUrso/GEEScripts:ProcessingModules/LC08_RockMask');

/*******************************************************
 * AOI
 *******************************************************/
var aoi_fc   = ee.FeatureCollection(AOI_ASSET);
var aoi_geom = aoi_fc.geometry();
Map.addLayer(aoi_geom, {color: '1f1e25'}, 'AOI');

/*******************************************************
 * METADATA TO KEEP
 *******************************************************/
var PROPS_TO_COPY = [
  'LANDSAT_PRODUCT_ID',
  'LANDSAT_SCENE_ID',
  'COLLECTION_CATEGORY',
  'SENSOR_ID',
  'SPACECRAFT_ID',
  'WRS_PATH',
  'WRS_ROW',
  'system:time_start',
  'CLOUD_COVER',
  'SUN_ELEVATION'
];

/*******************************************************
 * LOAD SCENES FROM LIST OF IDs
 *******************************************************/
var srcImages = SCENE_IDS.map(function(id) {
  var mission = id.slice(0, 4);   // 'LC08' or 'LC09'
  var colId   = 'LANDSAT/' + mission + '/C02/T1_TOA/' + id;

  return ee.Image(colId)
    .set('LANDSAT_PRODUCT_ID', id)
    .set('SENSOR_ID', mission);
});

var srcIC = ee.ImageCollection(srcImages);
print('Loaded scenes:', srcIC);
print('Source collection size:', srcIC.size());

/*******************************************************
 * add NDWI_ICE
 *******************************************************/
// NDWI_ICE = (B2 - B4) / (B2 + B4)
function addNdwiIce(img) {
  var ndwiIce = ee.Image(img)
    .normalizedDifference(['B2', 'B4'])
    .rename('NDWI_ICE')
    .toFloat();

  return ee.Image(img)
    .addBands(ndwiIce);
}

/*******************************************************
 * apply masks + clip + add NDWI_ICE
 *******************************************************/
function applyMasksClipAndNdwi(img) {
  // Cloud mask
  img = cloud_mask.maskClouds(img);

  // Rock mask
  img = rock_mask.rock_mask(img);

  // Clip to AOI
  img = img.clip(aoi_geom);

  // Add NDWI_ICE band
  img = addNdwiIce(img);

  // Copy over properties (and keep everything else)
  img = img.copyProperties(img, PROPS_TO_COPY);

  return img;
}

/*******************************************************
 * PROCESS COLLECTION
 *******************************************************/
var procIC = srcIC.map(applyMasksClipAndNdwi);
print('Processed collection size:', procIC.size());

/*******************************************************
 * OPTIONAL MAP PREVIEW
 *******************************************************/
if (ADD_RGB_LAYERS) {
  var list = procIC.toList(procIC.size());
  var pids = procIC.aggregate_array('LANDSAT_PRODUCT_ID');

  pids.evaluate(function(pidArray) {
    pidArray.forEach(function(pid, i) {
      var img = ee.Image(list.get(i));
      Map.addLayer(img, VIZ_RGB, 'RGB ' + pid);
      Map.addLayer(img.select('NDWI_ICE'), VIZ_NDWI_ICE, 'NDWI_ICE ' + pid, false);
    });
  });
}

/*******************************************************
 * EXPORTS: queue one task per image
 *******************************************************/
var procList   = procIC.toList(procIC.size());
var productIds = procIC.aggregate_array('LANDSAT_PRODUCT_ID').getInfo();

print('Product IDs for export:', productIds);

for (var i = 0; i < productIds.length; i++) {
  var pid  = productIds[i] || ('scene_' + i);
  var name = pid.replace(/\s+/g, '_');

  var img = ee.Image(procList.get(i));

  if (EXPORT_MODE === 'DRIVE') {
    //cast to Float32 to avoid mixedtype export errors
    img = img.toFloat();

    Export.image.toDrive({
      image: img,
      description: NAME_PREFIX + name,
      folder: DRIVE_FOLDER,
      fileNamePrefix: NAME_PREFIX + name,
      region: aoi_geom,
      scale: EXPORT_SCALE_M,
      maxPixels: MAX_PIXELS
     
      // crs: 'EPSG:3413',
      // fileFormat: 'GeoTIFF'
    });

  } else if (EXPORT_MODE === 'ASSET') {
    // (DONT FORGET TO make sure DEST_ASSET_PREFIX ends with '/')
    var assetId = DEST_ASSET_PREFIX + (NAME_PREFIX + name);

    Export.image.toAsset({
      image: img,                 
      description: NAME_PREFIX + name,
      assetId: assetId,
      region: aoi_geom,
      scale: EXPORT_SCALE_M,
      maxPixels: MAX_PIXELS
    });

  } else {
    print('EXPORT_MODE must be "DRIVE" or "ASSET". Current:', EXPORT_MODE);
  }
}
