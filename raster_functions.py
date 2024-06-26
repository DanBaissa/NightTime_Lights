import ssl
import sys
import shutil
import threading
from tkinter import Tk, Button, Label, Entry, StringVar, filedialog, messagebox, scrolledtext, Checkbutton, IntVar, ttk
from urllib.request import urlopen, Request
from datetime import datetime
from countries import load_countries, get_bounding_box
from nasa_cmr_api import search_nasa_cmr
from process_files import process_h5_files
from rasterio.mask import mask
from shapely.geometry import mapping
import glob
from astropy.stats import SigmaClip
from tkinter import filedialog, messagebox, Label, Button, Entry, StringVar
import geopandas as gpd
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.io import MemoryFile
import rasterstats
import numpy as np
import os
import pandas as pd
from shapely.geometry import Polygon, Point
from geopandas.tools import sjoin
import matplotlib.pyplot as plt
from raster_cropper import crop_raster_with_shapefile, crop_rasters
from downloader import download, start_download_thread
from raster_cropper import crop_raster_with_shapefile, crop_rasters


def crop_processed_files(self):
    input_folder = self.destination_folder.get()
    country_name = self.country_var.get().strip()

    crop_rasters(input_folder, country_name)

def perform_stacking(self):
        input_folder = self.input_folder_var.get()
        output_folder = self.output_folder_var.get()
        sigma_value = float(self.sigma_var.get())

        # Specify the output filenames
        output_filename_mean = 'output_mean.tif'
        output_filename_sigma = 'output_sigma_clipped.tif'
        output_plot = 'output_plot.pdf'

        # Get a list of all GeoTIFF files in the folder
        file_list = glob.glob(os.path.join(input_folder, '*.tif'))

        # Initialize a list to store the data arrays from each GeoTIFF file
        data_arrays = []

        for file in file_list:
            with rasterio.open(file) as src:
                data = src.read(1).astype('float32')  # Assuming data is in the first band
                data[data > 10000] = np.nan  # Example threshold, adjust as needed
                data_arrays.append(data)

        if not data_arrays:
            messagebox.showerror("Error", "No TIFF files found in the input folder.")
            return

        # Stack the arrays along a new dimension (axis 0)
        stacked_arrays = np.stack(data_arrays, axis=0)

        # Calculate the mean along axis 0 while ignoring NaN values
        mean_array = np.nanmean(stacked_arrays, axis=0)

        # Perform sigma clipping
        sigma_clip = SigmaClip(sigma=sigma_value)
        clipped_arrays = sigma_clip(stacked_arrays, axis=0)

        # Calculate the mean of the clipped arrays
        mean_array_clipped = np.nanmean(clipped_arrays, axis=0)

        # Get the metadata from the first GeoTIFF file
        with rasterio.open(file_list[0]) as src:
            meta = src.meta

        # Update the metadata
        meta.update(count=1)

        # Write the mean array to a new GeoTIFF file
        with rasterio.open(os.path.join(output_folder, output_filename_mean), 'w', **meta) as dst:
            dst.write(mean_array, 1)

        # Write the sigma clipped mean array to a new GeoTIFF file
        with rasterio.open(os.path.join(output_folder, output_filename_sigma), 'w', **meta) as dst:
            dst.write(mean_array_clipped, 1)

        # Plot the results
        fig, axs = plt.subplots(1, 2, figsize=(10, 5))
        axs[0].imshow(np.log1p(mean_array), cmap='turbo')
        axs[0].set_title('Log of Mean Stack')
        axs[1].imshow(np.log1p(mean_array_clipped), cmap='turbo')
        axs[1].set_title('Log of Sigma Clipped')
        plt.savefig(os.path.join(output_folder, output_plot))

        messagebox.showinfo("Success", "Stacking process completed.")

def process_raster_data(self, selected_country, raster_paths, csv_path, output_dir):
        # Filter the shapefile for the selected country
        gdf = self.shapefile[self.shapefile['COUNTRY'] == selected_country]

        # Reproject the country shape to match the CRS of the raster data (EPSG:3857)
        gdf = gdf.to_crs(epsg=3857)

        if gdf.empty:
            messagebox.showerror("Error", f"No shape data found for the selected country: {selected_country}")
            return

        # Crop rasters to the country boundary
        for raster_path in raster_paths:
            cropped_raster_path = os.path.join(output_dir, os.path.basename(raster_path))
            crop_raster_with_shapefile(raster_path, gdf, cropped_raster_path)

        # Create grid
        xmin, ymin, xmax, ymax = gdf.total_bounds
        width = height = 10_000  # 10km grid size
        rows = int(np.ceil((ymax - ymin) / height))
        cols = int(np.ceil((xmax - xmin) / width))

        polygons = []
        for x in range(cols):
            for y in range(rows):
                polygons.append(Polygon([(x * width + xmin, y * height + ymin),
                                         ((x + 1) * width + xmin, y * height + ymin),
                                         ((x + 1) * width + xmin, (y + 1) * height + ymin),
                                         (x * width + xmin, (y + 1) * height + ymin)]))

        grid = gpd.GeoDataFrame({'geometry': polygons}, crs='EPSG:3857')
        intersection = gpd.overlay(gdf, grid, how='intersection')

        # Process rasters
        for raster_path in raster_paths:
            with rasterio.open(raster_path) as src:
                nodata_value = src.nodata if src.nodata else -32768
                raster_data = src.read(1)
                raster_data = np.where(raster_data == nodata_value, np.nan, raster_data)  # replace nodata with NaN

                transform, width, height = calculate_default_transform(src.crs, intersection.crs, src.width, src.height,
                                                                       *src.bounds)
                kwargs = src.meta.copy()
                kwargs.update({
                    'crs': intersection.crs,
                    'transform': transform,
                    'width': width,
                    'height': height
                })

                with MemoryFile() as memfile:
                    with memfile.open(**kwargs) as mem_dst:
                        reproject(
                            source=rasterio.band(src, 1),
                            destination=rasterio.band(mem_dst, 1),
                            src_transform=src.transform,
                            src_crs=src.crs,
                            dst_transform=transform,
                            dst_crs=intersection.crs,
                            resampling=Resampling.nearest
                        )
                    reprojected_raster = memfile.open().read(1)

                zonal_stats = rasterstats.zonal_stats(intersection, reprojected_raster, affine=transform, stats='mean')
                raster_column_name = os.path.splitext(os.path.basename(raster_path))[0]
                intersection[raster_column_name] = [stat['mean'] for stat in zonal_stats]

        # Process CSV only if a path is provided
        if csv_path:
            df = pd.read_csv(csv_path)
            gdf_points = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude))
            gdf_points.set_crs("EPSG:4326", inplace=True)
            gdf_points = gdf_points.to_crs(epsg=3857)

            joined = sjoin(gdf_points, intersection, how='inner')
            for column in df.columns.difference(['longitude', 'latitude']):
                mean_obs = joined.groupby('index_right')[column].mean()
                intersection[column] = mean_obs

        intersection.to_file(os.path.join(output_dir, 'intersection.shp'))

        # Visualization
        for column in [col for col in intersection.columns if col not in ['geometry']]:
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
            intersection.plot(column=column, ax=ax, legend=True, edgecolor='black')
            plt.title(f'{column}')
            plt.savefig(os.path.join(output_dir, f'{column}.png'))