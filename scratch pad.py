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
from raster_functions import crop_processed_files, perform_stacking

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



    def browse_destination_folder(self):
        self.destination_folder.set(filedialog.askdirectory())

    def browse_input_folder(self):
        self.input_folder_var.set(filedialog.askdirectory())

    def browse_output_folder(self):
        self.output_folder_var.set(filedialog.askdirectory())







if __name__ == "__main__":
    shapefile_path = 'Data/Black_Marble_IDs/Black_Marble_World_tiles.shp'
    boundary_shapefile_path = 'Data/Black_Marble_IDs/BlackMarbleTiles/BlackMarbleTiles.shp'
    shapefile = load_countries(shapefile_path)
    app = URLSearcher()
    app.mainloop()
