import os
import rasterio
from rasterio.mask import mask
from shapely.geometry import mapping
import geopandas as gpd
import tempfile

def crop_raster_with_shapefile(raster_path, country_shapefile_path, selected_country, output_path):
    # Load the shapefile
    country_shape = gpd.read_file(country_shapefile_path)
    country_shape = country_shape[country_shape['COUNTRY'] == selected_country]

    with rasterio.open(raster_path) as src:
        # Ensure the CRS of the raster and shapefile match
        country_shape = country_shape.to_crs(src.crs)

        # Create a temporary raster file with the adjusted data
        with tempfile.NamedTemporaryFile(delete=False, suffix='.tif') as temp_file:
            with rasterio.open(temp_file.name, 'w', **src.meta) as temp_dst:
                temp_dst.write(src.read(1), 1)

            # Crop the adjusted raster with the shapefile
            with rasterio.open(temp_file.name) as temp_src:
                out_image, out_transform = mask(temp_src, [mapping(geom) for geom in country_shape.geometry], crop=True, invert=False)
                out_meta = temp_src.meta.copy()

        # Update the metadata
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })

        # Create directory for output file if it doesn't exist
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # Write the cropped raster to a new file
        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(out_image[0], 1)  # Select the first band

        # Delete the temporary file
        os.remove(temp_file.name)
