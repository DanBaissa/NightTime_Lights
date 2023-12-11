import os
import ssl
import sys
import shutil
import threading
from tkinter import Tk, Button, Label, Entry, StringVar, filedialog, messagebox, scrolledtext, ttk
from urllib.request import urlopen, Request
from datetime import datetime
from countries import load_countries, get_bounding_box
from nasa_cmr_api import search_nasa_cmr
from process_files import process_h5_files
import geopandas as gpd

USERAGENT = 'tis/download.py_1.0--' + sys.version.replace('\n','').replace('\r','')


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

        Button(self, text="Download", command=self.download).grid(row=6, column=0, columnspan=2)

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

        threading.Thread(target=self.start_download_thread, args=(urls, headers, CTX)).start()

    def start_download_thread(self, urls, headers, ctx):
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

        messagebox.showinfo('Download complete', 'Download of all files completed.')

        process_h5_files(self.country_var.get().strip(),
                         boundary_shapefile_path,
                         self.country_var.get().strip(),
                         self.destination_folder.get(),
                         self.destination_folder.get())

    def browse_destination_folder(self):
        self.destination_folder.set(filedialog.askdirectory())

if __name__ == "__main__":
    shapefile_path = 'Data/Black_Marble_IDs/Black_Marble_World_tiles.shp'
    boundary_shapefile_path = 'Data/Black_Marble_IDs/BlackMarbleTiles/BlackMarbleTiles.shp'
    shapefile = load_countries(shapefile_path)
    app = URLSearcher()
    app.mainloop()
