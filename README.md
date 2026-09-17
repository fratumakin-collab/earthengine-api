/***************************************************************
 ANALISIS NDBI NEGERI SOUHOKU, KECAMATAN AMAHAI
 TAHUN 2024 (LANDSAT 8 OLI/TIRS)

 Resolusi spasial: 30 meter
 Rumus: NDBI = (Float B7 - Float B5) / (Float B7 + Float B5)
 (Band 7 = SWIR 2, Band 5 = NIR)
***************************************************************/

// ============================================================
// 1. KOORDINAT BATAS NEGERI SOUHOKU
// ============================================================
var souhokuCoords = [
  [128.925557, -3.337351], [128.925557, -3.338063], [128.925557, -3.338661],
  [128.925899, -3.338661], [128.925899, -3.341303], [128.925566, -3.341303],
  [128.925563, -3.342269], [128.925170, -3.342269], [128.925170, -3.342718],
  [128.926615, -3.343118], [128.932951, -3.344873], [128.934763, -3.345375],
  [128.944947, -3.346915], [128.938815, -3.349618], [128.938316, -3.349838],
  [128.937838, -3.351012], [128.937218, -3.352109], [128.935792, -3.354631],
  [128.934597, -3.356744], [128.933293, -3.356007], [128.932993, -3.355837],
  [128.931595, -3.355019], [128.930134, -3.353844], [128.928610, -3.352669],
  [128.926648, -3.350655], [128.923349, -3.347067], [128.922640, -3.346228],
  [128.922640, -3.345640], [128.922598, -3.345242], [128.921972, -3.344507],
  [128.920949, -3.343458], [128.920281, -3.343185], [128.918903, -3.342073],
  [128.917526, -3.340835], [128.915960, -3.339199], [128.914917, -3.338024],
  [128.913915, -3.336723], [128.913059, -3.335632], [128.912600, -3.334835],
  [128.912495, -3.334205], [128.912182, -3.333492], [128.911890, -3.332988],
  [128.911431, -3.332653], [128.911264, -3.332065], [128.911093, -3.331678],
  [128.910867, -3.331163], [128.910742, -3.330450], [128.910826, -3.330072],
  [128.910575, -3.329611], [128.909970, -3.328834], [128.910033, -3.328373],
  [128.910220, -3.327974], [128.910575, -3.327869], [128.911556, -3.327911],
  [128.912454, -3.328163], [128.913540, -3.328814], [128.914291, -3.329380],
  [128.914604, -3.329926], [128.914557, -3.330385], [128.914562, -3.330618],
  [128.914583, -3.331226], [128.914615, -3.331646], [128.914792, -3.332045],
  [128.915084, -3.332233], [128.915710, -3.332359], [128.916190, -3.332590],
  [128.916399, -3.332779], [128.916587, -3.333052], [128.916754, -3.333304],
  [128.916963, -3.333891], [128.917297, -3.334416], [128.917317, -3.334793],
  [128.917255, -3.335150], [128.918465, -3.336640], [128.919885, -3.338087],
  [128.920553, -3.337961], [128.920857, -3.337989], [128.921492, -3.338045],
  [128.921952, -3.337794], [128.922578, -3.337836], [128.923308, -3.337920],
  [128.924237, -3.337808], [128.924841, -3.337576], [128.925064, -3.337430],
  [128.925259, -3.337301], [128.925557, -3.337351]
];

var roi = ee.Geometry.Polygon([souhokuCoords], null, false);
Map.centerObject(roi, 14);

var batasROI = ee.Image().byte().paint({featureCollection: roi, color: 1, width: 2});
Map.addLayer(batasROI, {palette: ['FFFFFF']}, 'Batas Negeri Souhoku', true);

// ============================================================
// 2. FUNGSI MASKING AWAN LANDSAT 8
// ============================================================
function maskLandsat8(image) {
  var qaPixel = image.select('QA_PIXEL');
  var maskCloud = qaPixel.bitwiseAnd(1 << 3).eq(0);
  var maskCloudShadow = qaPixel.bitwiseAnd(1 << 4).eq(0);
  var maskCirrus = qaPixel.bitwiseAnd(1 << 2).eq(0);
  var maskFill = qaPixel.bitwiseAnd(1 << 0).eq(0);
  var cloudMask = maskFill.and(maskCloud).and(maskCloudShadow).and(maskCirrus);

  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  return image.addBands(opticalBands, null, true).updateMask(cloudMask);
}

// ============================================================
// 3. PEMANGGILAN CITRA TAHUN 2024
// ============================================================
var composite2024 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filterDate('2024-01-01', '2024-12-31')
  .map(maskLandsat8)
  .median()
  .clip(roi);

// Visualisasi Warna Alami Landsat 8 (RGB: B4, B3, B2)
Map.addLayer(composite2024, {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 0.0, max: 0.3, gamma: 1.2}, 'Landsat 8 RGB 2024', false);

// ============================================================
// 4. MENGHITUNG NDBI = (Float B7 - Float B5) / (Float B7 + Float B5)
// ============================================================
var band7 = composite2024.select('SR_B7').toFloat();
var band5 = composite2024.select('SR_B5').toFloat();
var pembilang = band7.subtract(band5);
var penyebut = band7.add(band5);
var maskPenyebut = penyebut.abs().gt(0.000001);

var ndbi2024 = pembilang.divide(penyebut).updateMask(maskPenyebut).rename('NDBI').clip(roi);

// ============================================================
// 5. KLASIFIKASI KEPADATAN BANGUNAN TAHUN 2024
// ============================================================
var rentangValid = ndbi2024.gte(-1).and(ndbi2024.lte(0.3));

var kelasNDBI = ee.Image(0)
  .where(ndbi2024.gte(-1).and(ndbi2024.lte(0)), 1)
  .where(ndbi2024.gt(0).and(ndbi2024.lte(0.1)), 2)
  .where(ndbi2024.gt(0.1).and(ndbi2024.lte(0.2)), 3)
  .where(ndbi2024.gt(0.2).and(ndbi2024.lte(0.3)), 4)
  .updateMask(rentangValid)
  .rename('Kelas_NDBI')
  .toByte()
  .clip(roi);

var warnaKelas = ['0066FF', 'FFFF00', 'FF9900', 'FF0000'];
Map.addLayer(kelasNDBI, {min: 1, max: 4, palette: warnaKelas}, 'Klasifikasi NDBI Tahun 2024', true);
Map.addLayer(batasROI, {palette: ['FFFFFF']}, 'Garis Batas Souhoku', true);

// ============================================================
// 6. PERHITUNGAN LUAS PENGGUNAAN LAHAN TAHUN 2024
// ============================================================
var idKelas = ee.List([1, 2, 3, 4]);
var rentangKelas = ee.List(['-1 sampai 0', '> 0 sampai 0,1', '> 0,1 sampai 0,2', '> 0,2 sampai 0,3']);
var namaKelas = ee.List([
  'Non-permukiman / Lahan tidak terbangun / Badan air atau vegetasi',
  'Permukiman jarang / Lahan terbangun kurang rapat',
  'Permukiman rapat / Lahan terbangun rapat',
  'Permukiman sangat rapat / Lahan terbangun sangat rapat'
]);
var warnaHex = ee.List(['#0066FF', '#FFFF00', '#FF9900', '#FF0000']);

var pixelAreaHa = ee.Image.pixelArea().divide(10000).rename('Luas_ha');
var totalLuasValid = ee.Number(pixelAreaHa.updateMask(rentangValid).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: roi, scale: 30, maxPixels: 1e13, tileScale: 4
}).get('Luas_ha'));

var tabelLuas = ee.FeatureCollection(idKelas.map(function(kelas) {
  kelas = ee.Number(kelas);
  var luas = ee.Number(pixelAreaHa.updateMask(kelasNDBI.eq(kelas)).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: roi, scale: 30, maxPixels: 1e13, tileScale: 4
  }).get('Luas_ha', 0));
  var persentase = ee.Number(ee.Algorithms.If(totalLuasValid.gt(0), luas.divide(totalLuasValid).multiply(100), 0));
  var idx = kelas.subtract(1);
  return ee.Feature(null, {
    Kelas: kelas,
    Rentang_NDBI: rentangKelas.get(idx),
    Klasifikasi: namaKelas.get(idx),
    Warna: warnaHex.get(idx),
    Luas_ha: luas,
    Persentase: persentase,
    Tahun: 2024,
    Sensor: 'Landsat 8 OLI'
  });
}));

print('TABEL LUAS PENGGUNAAN LAHAN TAHUN 2024:', tabelLuas.sort('Kelas'));

// Grafik Luas
print(ui.Chart.feature.byFeature(tabelLuas.sort('Kelas'), 'Rentang_NDBI', ['Luas_ha'])
  .setChartType('ColumnChart').setOptions({
    title: 'Luas Penggunaan Lahan Berdasarkan NDBI Tahun 2024 (Landsat 8)',
    hAxis: {title: 'Rentang NDBI'}, vAxis: {title: 'Luas (Ha)'}, colors: ['#E65100']
  }));

// ============================================================
// 7. PANEL LEGENDA DAN TABEL HASIL DI LAYAR PETA
// ============================================================
var panel = ui.Panel({style: {position: 'bottom-left', width: '510px', padding: '10px', backgroundColor: 'white'}});
panel.add(ui.Label({value: 'NDBI NEGERI SOUHOKU TAHUN 2024', style: {fontWeight: 'bold', fontSize: '16px'}}));
panel.add(ui.Label({value: 'Sensor: Landsat 8 OLI/TIRS | Resolusi: 30 m', style: {fontSize: '11px', color: '#555'}}));
panel.add(ui.Label({value: 'Rumus: NDBI = (Float B7 - Float B5) / (Float B7 + Float B5)', style: {fontSize: '11px', color: '#555'}}));

function buatBaris(col, txt) {
  return ui.Panel({
    widgets: [
      ui.Label({style: {backgroundColor: col, padding: '8px', margin: '0 6px 4px 0', border: '1px solid #222'}}),
      ui.Label({value: txt, style: {fontSize: '11px', margin: '2px 0 0 0'}})
    ],
    layout: ui.Panel.Layout.flow('horizontal')
  });
}
panel.add(buatBaris('#0066FF', 'Kelas 1 (-1 s/d 0) : Non-permukiman / Air / Vegetasi'));
panel.add(buatBaris('#FFFF00', 'Kelas 2 (>0 s/d 0,1) : Permukiman Jarang'));
panel.add(buatBaris('#FF9900', 'Kelas 3 (>0,1 s/d 0,2) : Permukiman Rapat'));
panel.add(buatBaris('#FF0000', 'Kelas 4 (>0,2 s/d 0,3) : Permukiman Sangat Rapat'));

panel.add(ui.Label({value: 'Rincian Luas Penggunaan Lahan (2024):', style: {fontWeight: 'bold', fontSize: '13px', margin: '8px 0 4px 0'}}));

var header = ui.Panel({
  widgets: [
    ui.Label({value: 'Kelas', style: {fontWeight: 'bold', width: '45px'}}),
    ui.Label({value: 'Rentang', style: {fontWeight: 'bold', width: '110px'}}),
    ui.Label({value: 'Luas (Ha)', style: {fontWeight: 'bold', width: '100px'}}),
    ui.Label({value: 'Persentase', style: {fontWeight: 'bold', width: '90px'}})
  ], layout: ui.Panel.Layout.flow('horizontal')
});
panel.add(header);

tabelLuas.sort('Kelas').evaluate(function(res) {
  res.features.forEach(function(f) {
    var p = f.properties;
    panel.add(ui.Panel({
      widgets: [
        ui.Label({value: String(p.Kelas), style: {width: '45px'}}),
        ui.Label({value: p.Rentang_NDBI, style: {width: '110px'}}),
        ui.Label({value: Number(p.Luas_ha).toFixed(2), style: {width: '100px'}}),
        ui.Label({value: Number(p.Persentase).toFixed(2) + ' %', style: {width: '90px'}})
      ], layout: ui.Panel.Layout.flow('horizontal')
    }));
  });
  totalLuasValid.evaluate(function(tot) {
    panel.add(ui.Label({value: 'Total Luas: ' + Number(tot).toFixed(2) + ' Ha', style: {fontWeight: 'bold', margin: '6px 0 0 0'}}));
  });
});
Map.add(panel);

// ============================================================
// 8. EKSPOR KE GOOGLE DRIVE
// ============================================================
Export.table.toDrive({
  collection: tabelLuas.sort('Kelas'), description: 'Tabel_NDBI_Souhoku_2024',
  folder: 'NDBI_Souhoku_2024', fileNamePrefix: 'Tabel_NDBI_Souhoku_2024', fileFormat: 'CSV'
});

Export.image.toDrive({
  image: kelasNDBI, description: 'Klasifikasi_NDBI_Souhoku_2024',
  folder: 'NDBI_Souhoku_2024', fileNamePrefix: 'Klasifikasi_NDBI_Souhoku_2024',
  region: roi, scale: 30, crs: 'EPSG:32752', maxPixels: 1e13
});
