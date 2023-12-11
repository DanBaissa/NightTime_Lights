import geopandas as gpd

def load_countries(shapefile_path):
    gdf = gpd.read_file(shapefile_path)
    return gdf

def get_bounding_box(country, gdf):
    country_gdf = gdf[gdf['COUNTRY'] == country]
    return country_gdf.total_bounds
