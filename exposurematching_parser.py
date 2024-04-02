import tkinter as tk
from tkinter import filedialog
import rasterio
from skimage import exposure
import numpy as np
import os

def normalize(image):
    return (image - image.min()) / (image.max() - image.min())

def load_and_match_histograms():
    # File dialogs to select the multiple source, a single reference, and a directory for output files
    source_paths = filedialog.askopenfilenames(title="Select Source GeoTIFFs")  # Allows multiple selection
    reference_path = filedialog.askopenfilename(title="Select Reference GeoTIFF")
    output_directory = filedialog.askdirectory(title="Select Output Directory")

    if source_paths and reference_path and output_directory:
        # Load the reference GeoTIFF
        with rasterio.open(reference_path) as ref:
            image_reference = ref.read(1)  # Reading the first channel

        for source_path in source_paths:
            with rasterio.open(source_path) as src:
                image_source = src.read(1)  # Reading the first channel
                meta_source = src.meta

                # Apply histogram matching
                matched_image = exposure.match_histograms(image_source, image_reference)

                # Normalize the matched image
                normalized_image = normalize(matched_image)

                # Construct unique output paths
                base_name = os.path.splitext(os.path.basename(source_path))[0]
                matched_output_path = os.path.join(output_directory, f"{base_name}_matched.tif")
                normalized_output_path = os.path.join(output_directory, f"{base_name}_normalized.tif")

                # Save the matched and normalized images
                meta_source.update(count=1)  # Ensure metadata is updated for single channel
                with rasterio.open(matched_output_path, "w", **meta_source) as dest:
                    dest.write(matched_image, 1)

                with rasterio.open(normalized_output_path, "w", **meta_source) as dest:
                    dest.write(normalized_image, 1)

        print("Histogram matching and normalization completed for all selected files. Files saved.")


# Set up the GUI window
root = tk.Tk()
root.title("GeoTIFF Histogram Matcher")
root.geometry("300x150")

# Add a button to start the process
process_button = tk.Button(root, text="Load GeoTIFFs and Match Histograms", command=load_and_match_histograms)
process_button.pack(pady=20)

# Start the GUI event loop
root.mainloop()
