import os
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
import geopandas as gpd
import glob
import rasterio
from rasterio.mask import mask
from shapely.geometry import mapping

USERAGENT = 'tis/download.py_1.0--' + sys.version.replace('\n', '').replace('\r', '')

# Load the shapefile
shapefile = gpd.read_file('Data/Black_Marble_IDs/Black_Marble_World_tiles.shp')
boundary_shapefile = gpd.read_file('Data/Black_Marble_IDs/BlackMarbleTiles/BlackMarbleTiles.shp')


class URLSearcher(Tk):
    def __init__(self):
        super().__init__()

        self.title("Black Marble Night Lights Downloader")

        self.start_date = StringVar()
        self.end_date = StringVar()
        self.country_var = StringVar()
        self.token = StringVar()
        self.destination_folder = StringVar()
        self.crop_files_var = IntVar()

        Label(self, text="Start Date (e.g., 'June 1, 2022'):").grid(row=0, column=0)
        Entry(self, textvariable=self.start_date).grid(row=0, column=1)

        Label(self, text="End Date (e.g., 'June 2, 2022'):").grid(row=1, column=0)
        Entry(self, textvariable=self.end_date).grid(row=1, column=1)

        Label(self, text="Country:").grid(row=2, column=0)
        self.country_dropdown = ttk.Combobox(self, textvariable=self.country_var)
        self.country_dropdown['values'] = [str(i) for i in shapefile['COUNTRY'].unique().tolist()]
        self.country_dropdown.grid(row=2, column=1)

        Label(self, text='Destination folder:').grid(row=3, column=0)
        Entry(self, textvariable=self.destination_folder).grid(row=3, column=1)
        Button(self, text='Browse...', command=self.browse_destination_folder).grid(row=3, column=2)

        Label(self, text='Token:').grid(row=4, column=0)
        Entry(self, textvariable=self.token).grid(row=4, column=1)

        self.progress = ttk.Progressbar(self, orient='horizontal', mode='determinate')
        self.progress.grid(row=5, column=0, columnspan=3)

        Checkbutton(self, text="Crop Files after Processing", variable=self.crop_files_var).grid(row=6, column=0, columnspan=2)

        Button(self, text="Download", command=self.download).grid(row=7, column=0, columnspan=2)

    def browse_destination_folder(self):
        self.destination_folder.set(filedialog.askdirectory())

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
