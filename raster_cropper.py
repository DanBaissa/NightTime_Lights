# raster_cropper.py

import glob
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from shapely.geometry import mapping
import os
import tempfile

def crop_raster_with_shapefile(raster_path, country_shape, output_path):
    print('Opening raster file...')
    with rasterio.open(raster_path) as src:
        print('Ensuring CRS match...')
        country_shape = country_shape.to_crs(src.crs)

        print('Cropping raster file...')
        out_image, out_transform = mask(src, [mapping(geom) for geom in country_shape.geometry], crop=True, invert=False)
        out_meta = src.meta.copy()

        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    print('Writing cropped raster to new file...')
    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(out_image[0], 1)

    print('Finished cropping raster file.')

def crop_rasters(input_folder, output_folder, country_name):
    country_shapefile = 'Data/World_Countries/World_Countries_Generalized.shp'
    shapefile = gpd.read_file(country_shapefile)
    country = shapefile[shapefile['COUNTRY'] == country_name]

    raster_files = glob.glob(os.path.join(input_folder, '*.tif'))

    for raster_file in raster_files:
        output_file = os.path.join(output_folder, os.path.basename(raster_file))
        crop_raster_with_shapefile(raster_file, country, output_file)
