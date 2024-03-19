import sys
import netCDF4 as nc
import geopandas as gpd
import os
import pyproj
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import rioxarray
import xarray as xr
import re

def main(station_number):
    # Define the coordinate reference system (CRS) for the data
    osgb36 = "EPSG:27700"

    # Load the catchment shapefile
    sf = gpd.read_file("~/NC_intl/other_work/data/shapefiles/nrfa_all_catchments.shp", crs=osgb36)
    # Extract the station IDs from the shapefile
    sf_id = sf["STATION"].values

    # Read the eflag catchment eflag station IDs from Excel
    eflag_id = pd.read_excel("~/NC_intl/other_work/data/eFLaG_Station_Metadata.xlsx", sheet_name="River Flow")
    station_names = eflag_id["NAME"].values
    eflag_id = eflag_id["STATION"].values

    # Find matching indices of stations between shapefile and eflag IDs
    matching_indices = np.where(np.isin(sf_id, eflag_id))[0]

    # Define the directory containing the NetCDF files
    filepath = "/badc/deposited2021/chess-scape/data/rcp85_bias-corrected/01/daily/tas/"

    # Get a list of NetCDF file names and sort them
    n_dir = os.listdir(filepath)
    n_dir_= sorted(n_dir, key=lambda x: re.findall(r'\d{8}', x))

    # Calculate the station index based on the provided station number
    station_index = station_number - 1  # Adjust for zero-based indexing
    print(station_index)
    # Select the matching station's geometry
    sf1 = sf.iloc[matching_indices]
    sf1 = sf1.to_crs(osgb36)
    selected_polygon = sf1.iloc[station_index]
    station_value = selected_polygon['STATION']
    selected_gdf = gpd.GeoDataFrame({"geometry": [selected_polygon.geometry]}, crs=sf.crs)
    selected_gdf = selected_gdf.to_crs(osgb36)

    # Create an empty list to store spatial mean arrays
    spatial_means_list = []

    # Loop through the NetCDF files
    for ii in range(len(n_dir)):
        filename = os.path.join(filepath, n_dir[ii])
        nc_data = nc.Dataset(filename)
        tasv = nc_data.variables["tas"][:]
        latv = nc_data.variables["y"][:]
        lonv = nc_data.variables["x"][:]
        
        # Create an empty list to store spatial means for each time step
        sma = []
        for t_index in range(len(tasv)):
            selected_tas = tasv[t_index, :, :]

            tasv_xr = xr.DataArray(selected_tas, coords=[("lat", latv), ("lon", lonv)])
            tasv_xr.rio.write_crs(osgb36, inplace=True)
            tasv_xr.rio.set_spatial_dims(x_dim="lon", y_dim="lat")

            # Clip data to the selected geometry and calculate spatial mean
            clipped_tas_xr = tasv_xr.rio.clip(selected_gdf.geometry)
            spatial_mean = np.array(clipped_tas_xr.mean(dim=("lat", "lon")))
            sma.append(spatial_mean - 273.15)

        # Append spatial means for this file to the list
        spatial_means_list.append(sma)
        rt = np.array(spatial_means_list)
        # print(spatial_means_list)

    # Reshape and assign spatial means for this station ID to the final array
    final_reshaped_array = rt.flatten()
    final_reshaped_array = np.round(final_reshaped_array, decimals=4)
    fnw = f"/gws/nopw/j04/nzplus/3A/srijha/data/chess1k/tas01/temp_daily{station_value}.csv"
    np.savetxt(fnw, final_reshaped_array, fmt='%.4f')

if __name__ == "__main__":
    # Check if the station number is provided as an argument
    if len(sys.argv) < 2:
        print("Usage: python script.py <station_number>")
        sys.exit(1)
    
    # Get the station number from command-line arguments
    station_number = int(sys.argv[1])
    
    # Call the main function with the station number
    main(station_number)
