import tkinter as tk
from tkinter import filedialog
import rasterio
from skimage import exposure
import numpy as np

def normalize(image):
    return (image - image.min()) / (image.max() - image.min())

def load_and_match_histograms():
    # File dialogs to select the source, reference, and output file paths
    source_path = filedialog.askopenfilename(title="Select Source GeoTIFF")
    reference_path = filedialog.askopenfilename(title="Select Reference GeoTIFF")
    output_path = filedialog.asksaveasfilename(defaultextension=".tif", title="Save Matched GeoTIFF")

    if source_path and reference_path and output_path:
        # Load the GeoTIFFs
        with rasterio.open(source_path) as src:
            image_source = src.read(1)  # Reading the first channel
            meta_source = src.meta

        with rasterio.open(reference_path) as ref:
            image_reference = ref.read(1)  # Reading the first channel

        # Apply histogram matching
        matched_image = exposure.match_histograms(image_source, image_reference)

        # Normalize the matched image
        normalized_image = normalize(matched_image)

        # Save the matched image
        meta_source.update(count=1)  # Ensure metadata is updated for single channel
        with rasterio.open(output_path, "w", **meta_source) as dest:
            dest.write(matched_image, 1)

        # Save the normalized image
        normalized_output_path = output_path.replace(".tif", "_normalized.tif")
        with rasterio.open(normalized_output_path, "w", **meta_source) as dest:
            dest.write(normalized_image, 1)

        print("Histogram matching and normalization completed. Files saved.")

# Set up the GUI window
root = tk.Tk()
root.title("GeoTIFF Histogram Matcher")
root.geometry("300x150")

# Add a button to start the process
process_button = tk.Button(root, text="Load GeoTIFFs and Match Histograms", command=load_and_match_histograms)
process_button.pack(pady=20)

# Start the GUI event loop
root.mainloop()
