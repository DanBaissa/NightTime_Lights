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

USERAGENT = 'tis/download.py_1.0--' + sys.version.replace('\n', '').replace('\r', '')

# Load the shapefile
shapefile = gpd.read_file('Data/Black_Marble_IDs/Black_Marble_World_tiles.shp')
boundary_shapefile = gpd.read_file('Data/Black_Marble_IDs/BlackMarbleTiles/BlackMarbleTiles.shp')


class URLSearcher(Tk):
    def __init__(self):
        super().__init__()

        self.title("Black Marble Night Lights Downloader")

        # Create a label with "Step 1: Downloading" and use grid
        step1_label = Label(self, text="Step 1: Downloading", font=("Arial", 14))
        step1_label.grid(row=0, column=0, columnspan=3, pady=(10, 0),
                         sticky='ew')  # Adjust row and columnspan as needed

        # Create a horizontal separator line and use grid
        separator = ttk.Separator(self, orient='horizontal')
        separator.grid(row=1, column=0, columnspan=3, sticky='ew', padx=5,
                       pady=(0, 10))  # Adjust row and columnspan as needed

        self.start_date = StringVar()
        self.end_date = StringVar()
        self.country_var = StringVar()
        self.token = StringVar()
        self.destination_folder = StringVar()
        self.crop_files_var = IntVar()

        Label(self, text="Start Date (e.g., 'June 1, 2022'):").grid(row=2, column=0)
        Entry(self, textvariable=self.start_date).grid(row=2, column=1)

        Label(self, text="End Date (e.g., 'June 2, 2022'):").grid(row=3, column=0)
        Entry(self, textvariable=self.end_date).grid(row=3, column=1)

        Label(self, text="Country:").grid(row=4, column=0)
        self.country_dropdown = ttk.Combobox(self, textvariable=self.country_var)
        self.country_dropdown['values'] = [str(i) for i in shapefile['COUNTRY'].unique().tolist()]
        self.country_dropdown.grid(row=4, column=1)

        Label(self, text='Destination folder:').grid(row=5, column=0)
        Entry(self, textvariable=self.destination_folder).grid(row=5, column=1)
        Button(self, text='Browse...', command=self.browse_destination_folder).grid(row=5, column=2)

        Label(self, text='Token:').grid(row=6, column=0)
        Entry(self, textvariable=self.token).grid(row=6, column=1)

        self.progress = ttk.Progressbar(self, orient='horizontal', mode='determinate')
        self.progress.grid(row=7, column=0, columnspan=3)

        Checkbutton(self, text="Crop Files after Processing", variable=self.crop_files_var).grid(row=8, column=0, columnspan=2)

        Button(self, text="Download", command=self.download).grid(row=9, column=0, columnspan=2)

        # Create a horizontal separator line and use grid
        separator = ttk.Separator(self, orient='horizontal')
        separator.grid(row=10, column=0, columnspan=3, sticky='ew', padx=0,
                       pady=(15, 0))  # Adjust row and columnspan as needed

        # Step 2: Stacking Data
        step2_label = Label(self, text="Step 2: Stacking Data", font=("Arial", 14))
        step2_label.grid(row=12, column=0, columnspan=3, pady=(10, 0), sticky='ew')

        separator2 = ttk.Separator(self, orient='horizontal')
        separator2.grid(row=13, column=0, columnspan=3, sticky='ew', padx=5, pady=(0, 10))

        self.input_folder_var = StringVar()
        self.output_folder_var = StringVar()
        self.sigma_var = StringVar(value="3")  # Default sigma value

        Label(self, text="Input Folder:").grid(row=14, column=0)
        Entry(self, textvariable=self.input_folder_var).grid(row=14, column=1)
        Button(self, text="Browse...", command=self.browse_input_folder).grid(row=14, column=2)

        Label(self, text="Output Folder:").grid(row=15, column=0)
        Entry(self, textvariable=self.output_folder_var).grid(row=15, column=1)
        Button(self, text="Browse...", command=self.browse_output_folder).grid(row=15, column=2)

        Label(self, text="Sigma Clipping:").grid(row=16, column=0)
        Entry(self, textvariable=self.sigma_var).grid(row=16, column=1)

        Button(self, text="Perform Stacking", command=self.perform_stacking).grid(row=17, column=0, columnspan=3)

        # Step 3: Raster Merger (New Functionality)
        step3_label = Label(self, text="Step 3: Raster Merger", font=("Arial", 14))
        step3_label.grid(row=18, column=0, columnspan=3, pady=(10, 0), sticky='ew')

        # Variables for file paths
        self.load_country_shapefile()
        self.country_var = StringVar()
        self.raster_paths_var = StringVar()
        self.csv_path_var = StringVar()
        self.output_dir_var = StringVar()

        # Country selection
        Label(self, text="Country:").grid(row=19, column=0)
        self.country_dropdown = ttk.Combobox(self, textvariable=self.country_var)
        self.country_dropdown['values'] = [str(country) for country in self.countries]
        self.country_dropdown.grid(row=19, column=1)

        # Raster selection
        Label(self, text="Rasters:").grid(row=20, column=0)
        Entry(self, textvariable=self.raster_paths_var).grid(row=20, column=1)
        Button(self, text="Browse...", command=self.browse_rasters).grid(row=20, column=2)

        # CSV file selection (optional)
        Label(self, text="CSV File (Optional):").grid(row=21, column=0)
        Entry(self, textvariable=self.csv_path_var).grid(row=21, column=1)
        Button(self, text="Browse...", command=self.browse_csv).grid(row=21, column=2)

        # Output directory selection
        Label(self, text="Output Directory:").grid(row=22, column=0)
        Entry(self, textvariable=self.output_dir_var).grid(row=22, column=1)
        Button(self, text="Browse...", command=self.browse_output_dir).grid(row=22, column=2)

        # Run button for Step 3
        self.run_step3_button = Button(self, text="Run Step 3", command=self.run_step3)
        self.run_step3_button.grid(row=23, column=0, columnspan=3)

    # New methods for file/folder browsing...
    def browse_shapefile(self):
        self.shapefile_path_var.set(
            filedialog.askopenfilename(filetypes=(("Shapefile", "*.shp"), ("All files", "*.*"))))

    def browse_rasters(self):
        self.raster_paths_var.set(filedialog.askopenfilenames(filetypes=(("GeoTIFF", "*.tif"), ("All files", "*.*"))))

    def browse_csv(self):
        self.csv_path_var.set(filedialog.askopenfilename(filetypes=(("CSV Files", "*.csv"), ("All files", "*.*"))))

    def browse_output_dir(self):
        self.output_dir_var.set(filedialog.askdirectory())

    # Method for running Step 3 functionality

    def browse_rasters(self):
        file_paths = filedialog.askopenfilenames(filetypes=[("GeoTIFF", "*.tif"), ("All files", "*.*")])
        if file_paths:  # Check if any file is selected
            self.raster_paths_var.set(','.join(file_paths))  # Join the paths with a comma

    def load_country_shapefile(self):
        country_shapefile = 'Data/World_Countries/World_Countries_Generalized.shp'
        self.shapefile = gpd.read_file(country_shapefile)
        self.countries = self.shapefile['COUNTRY'].unique()

    def run_step3(self):
        selected_country = self.country_var.get().strip()
        raster_paths = self.raster_paths_var.get().split(',')  # Split the paths by comma
        csv_path = self.csv_path_var.get()
        output_dir = self.output_dir_var.get()

        if not selected_country or not raster_paths or not output_dir:
            messagebox.showwarning("Warning", "Please select all required inputs and output directory.")
            return

        country_shape = self.shapefile[self.shapefile['COUNTRY'] == selected_country]

        try:
            self.process_raster_data(selected_country, raster_paths, csv_path, output_dir)
            messagebox.showinfo("Success", "Step 3 processing completed successfully.")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {e}")

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

    def browse_destination_folder(self):
        self.destination_folder.set(filedialog.askdirectory())

    def browse_input_folder(self):
        self.input_folder_var.set(filedialog.askdirectory())

    def browse_output_folder(self):
        self.output_folder_var.set(filedialog.askdirectory())

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


    def download(self):
        start_date_str = self.start_date.get()
        end_date_str = self.end_date.get()
        selected_country = self.country_var.get().strip()
        selected_rows = shapefile[shapefile['COUNTRY'] == selected_country]
        tiles = selected_rows['TileID'].tolist()

        start_date = datetime.strptime(start_date_str, '%B %d, %Y')
        end_date = datetime.strptime(end_date_str, '%B %d, %Y')
        bounding_box = get_bounding_box(selected_country, shapefile)
        collection_id = "C1898025206-LAADS"  # Example collection ID

        urls = search_nasa_cmr(collection_id, start_date, end_date, bounding_box)

        messagebox.showinfo("", "Download Started.")

        self.progress['maximum'] = len(urls)
        self.progress['value'] = 0

        CTX = ssl.SSLContext()
        headers = {'user-agent': USERAGENT, 'Authorization': 'Bearer ' + self.token.get()}

        crop_files = self.crop_files_var.get() == 1
        threading.Thread(target=self.start_download_thread, args=(urls, headers, CTX, crop_files)).start()

    def start_download_thread(self, urls, headers, ctx, crop_files):
        for url in urls:
            url = url.strip()
            filename = url.split('/')[-1]
            dest = os.path.join(self.destination_folder.get(), filename)

            try:
                with urlopen(Request(url, headers=headers), context=ctx) as response:
                    with open(dest, 'wb') as out_file:
                        shutil.copyfileobj(response, out_file)
                self.progress['value'] += 1
            except Exception as e:
                messagebox.showerror('Download error', f'Failed to download {url} due to {str(e)}')
                break

        process_h5_files(self.country_var.get().strip(),
                         boundary_shapefile_path,
                         self.country_var.get().strip(),
                         self.destination_folder.get(),
                         self.destination_folder.get())

        if crop_files:
            self.crop_processed_files()

        messagebox.showinfo('Download complete', 'Download of all files completed.')

    def crop_processed_files(self):
        input_folder = self.destination_folder.get()
        country_name = self.country_var.get().strip()

        crop_rasters(input_folder, country_name)

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

def crop_rasters(input_folder, country_name):
    country_shapefile = 'Data/World_Countries/World_Countries_Generalized.shp'
    shapefile = gpd.read_file(country_shapefile)
    country = shapefile[shapefile['COUNTRY'] == country_name]

    raster_files = glob.glob(os.path.join(input_folder, '*.tif'))

    for raster_file in raster_files:
        # The output file will be the same as the input file, effectively overwriting it
        crop_raster_with_shapefile(raster_file, country, raster_file)

if __name__ == "__main__":
    shapefile_path = 'Data/Black_Marble_IDs/Black_Marble_World_tiles.shp'
    boundary_shapefile_path = 'Data/Black_Marble_IDs/BlackMarbleTiles/BlackMarbleTiles.shp'
    shapefile = load_countries(shapefile_path)
    app = URLSearcher()
    app.mainloop()
