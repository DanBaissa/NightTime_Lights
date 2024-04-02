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