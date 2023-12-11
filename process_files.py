import os
import re
import h5py
import rasterio
import contextlib
from rasterio.merge import merge
from rasterio.transform import from_bounds
import geopandas as gpd

# Load the shapefile
shapefile = gpd.read_file('Data/Black_Marble_IDs/Black_Marble_World_tiles.shp')
boundary_shapefile = gpd.read_file('Data/Black_Marble_IDs/BlackMarbleTiles/BlackMarbleTiles.shp')


def process_h5_files(country_shapefile_path, boundary_shapefile_path, selected_country, input_folder, output_folder):
    # Filter the shapefile to get the rows for the selected country
    selected_rows = shapefile[shapefile['COUNTRY'] == selected_country.strip()]

    # Get a list of all .h5 files in the selected directory
    files = [f for f in os.listdir(input_folder) if f.endswith('.h5')]

    # Create a dictionary where the keys are the days and the values are the lists of files for each day
    files_by_day = {}
    for file in files:
        day = re.search('\.A(\d+)\.', file).group(1)
        if day not in files_by_day:
            files_by_day[day] = []
        files_by_day[day].append(file)

    # Loop through the days and process the files for each day
    for day, files in files_by_day.items():
        # List to store the raster data
        rasters = []

        # Loop through the files and process them if their TileID is in selected_rows
        for file in files:
            tile_id = re.search('h\d{2}v\d{2}', file).group()
            shape = selected_rows[selected_rows['TileID'] == tile_id]

            if not shape.empty:
                with h5py.File(os.path.join(input_folder, file), "r") as f:
                    data_fields = f["/HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields"]
                    gap_filled_dnb_brdf_corrected_ntl = data_fields['Gap_Filled_DNB_BRDF-Corrected_NTL'][...]

                    # Get the bounds of the shape from the boundary shapefile
                    boundary_shape = boundary_shapefile[boundary_shapefile['TileID'] == tile_id]
                    left, bottom, right, top = boundary_shape.total_bounds

                    # Save the raster data to a temporary GeoTIFF file
                    temp_file = os.path.join(input_folder, f"temp_{tile_id}.tif")
                    with rasterio.open(
                            temp_file,
                            "w",
                            driver="GTiff",
                            height=gap_filled_dnb_brdf_corrected_ntl.shape[0],
                            width=gap_filled_dnb_brdf_corrected_ntl.shape[1],
                            count=1,
                            dtype=gap_filled_dnb_brdf_corrected_ntl.dtype,
                            crs=shape.crs.to_epsg(),  # Use the CRS of the shapefile
                            transform=from_bounds(left, bottom, right, top,
                                                  gap_filled_dnb_brdf_corrected_ntl.shape[1],
                                                  gap_filled_dnb_brdf_corrected_ntl.shape[0])
                    ) as new_dataset:
                        new_dataset.write(gap_filled_dnb_brdf_corrected_ntl, 1)

                    rasters.append(temp_file)

        # Open the rasters
        with contextlib.ExitStack() as stack:
            files_to_merge = [stack.enter_context(rasterio.open(raster)) for raster in rasters]

            # Merge the rasters
            merged, transform = merge(files_to_merge)

        # Remove the temporary files
        for raster in rasters:
            os.remove(raster)

        output_file = os.path.join(output_folder, f"{selected_country}_{day}.tif")
        with rasterio.open(
                output_file,
                "w",
                driver="GTiff",
                height=merged.shape[1],
                width=merged.shape[2],
                count=1,
                dtype=merged.dtype,
                crs=shape.crs.to_epsg(),
                transform=transform
        ) as dest:
            dest.write(merged)

            # Delete the .h5 files after processing
            for file in files:
                os.remove(os.path.join(input_folder, file))
