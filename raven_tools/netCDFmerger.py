import shutil

import numpy as np
from netCDF4 import Dataset
from numpy import exp


def compute_lat_radians(values):
    """Slow way of computing latitude in radians

    See https://jakevdp.github.io/PythonDataScienceHandbook/02.03-computation-on-arrays-ufuncs.html
    :param values:
    :return:
    """
    output = np.empty(len(values))
    for i in range(len(values)):
        output[i] = values[i] * np.pi / 180
    return output


if __name__ == '__main__':
    # Path to the clipped netCDF file to work with
    cdf_in_path = "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/TabsD_v2.0_swiss.lv95/out/TabsD_ch01r.swiss" \
                  ".lv95_200001010000_200012310000_clipped.nc"
    # Path to the netCDF file that will be appended with PET data
    cdf_out_path = "/media/mainman/Data/RAVEN/data/MeteoSwiss_gridded_products/TabsD_v2.0_swiss.lv95/out/TabsD_ch01r" \
                   ".swiss.lv95_200001010000_200012310000_pet.nc"
    # Copy the clipped file to a new file so as not to change the original file
    shutil.copyfile(cdf_in_path, cdf_out_path)

    # Read the clipped netCDF file into a netCDF Dataset
    cdf_dataset_in = Dataset(cdf_in_path, format="NETCDF4")
    # Create the output file in append mode
    cdf_dataset_out = Dataset(cdf_out_path, "r+", format="NETCDF4")

    # Create a new variable 'PET' in the output file, according to the 'TabsD' variable in the input file
    pet = cdf_dataset_out.createVariable('PET', np.float32, fill_value=-999.99, dimensions=('time', 'N', 'E'))
    pet.units = 'mm/d'
    pet.long_name = 'daily potential evapotranspiration (Hamon)'
    pet.setncatts({'grid_name': "ch01r.swiss.lv95",
                   'version': "v1.2",
                   'prod_date': "2022-10-01",
                   'coordinates': "lon lat",
                   'grid_mapping': u"swiss_lv95_coordinates"})
    # Get the values of 'TabsD' of the dataset as a nested array
    tabsd = cdf_dataset_out['TabsD'][:]
    latitude = cdf_dataset_out['lat'][:,1]
    pet_array = cdf_dataset_out['PET'][:]
    # Alternatively, Get the 'TabsD' variable as a netCDF4 Variable data type
    # tabsd = cdf_dataset_out.variables['TabsD']
    dl = np.empty((366,60,36))

    # for lx, lat in np.ndenumerate(latitude):
    #     latitude[lx] = lat * np.pi /180


    # for dx, tmean in np.ndenumerate(tabsd):
    #     print(f"Index: {dx}")
    #     print(f"TabsD: {tmean}")
    #     print(f"Lat(deg): {}")
    #     lrad[dx[0],dx[1],dx[2]] = dx[1] * np.pi /180
    #     print(f"Lat(rad): {lrad[dx[0],dx[1],dx[2]]}")
        # output_dl[dx[0],dx[1],dx[2]] = ((24 / np.pi * (np.arccos(-np.tan(0.409 * np.sin(2. * np.pi / 365. * (dx[0]+1) - 1.39)) * np.tan(dx[1] * np.pi /180)))) / 12) ** 2 * exp(day / 16)
        # print(output_dl[dx[0],dx[1],dx[2]])

    for dx, day in np.ndenumerate(tabsd):
        sol_dec = 0.409 * np.sin(2. * np.pi / 366. * (dx[0]+1) - 1.39)
        l_rad = latitude[dx[1]] * np.pi / 180
        s_angle = np.arccos(-np.tan(sol_dec) * np.tan(l_rad))
        dlh = 24 / np.pi * s_angle
        pet_value = (dlh / 12) ** 2 * exp(day / 16)
        # print(f"dx: {dx}")
        # print(f"lat(rad): {l_rad}")
        # print(f"sol_dec: {sol_dec}")
        # print(f"s_angle: {s_angle}")
        # print(f"Daylight hours: {dlh}")
        # print(f"PET(Hamon: {pet_value}")
        dl[dx[0],dx[1],dx[2]] = pet_value
        # print(f"PET(Hamon: {dl[dx[0],dx[1],dx[2]]}")
    cdf_dataset_out['PET'][:] = dl


    # Get the 'time' variable as a netCDF4 Variable data type
    time = cdf_dataset_out.variables['time']
    # From the 'time' netCDF4 Variable, get the values as a nested array
    time_array = time[:]
    tmean_array = tabsd[:]

    # Calculate the latitude in radians from degrees
    lat_radians = tmean_array * np.pi / 180

    # Get the values for day=1
    tabsd[1, :, :]

    # Another way of getting all the data from one netCDF file and writing it into a new one. However, it gave me an
    # error, since _FillValue cannot be written except when creating a variable with cdf_dataset_out as cdf_dataset_out:
    # copy global attributes all at once via dictionary # cdf_dataset_out.setncatts(cdf_dataset_in.__dict__) # copy
    # dimensions for name, dimension in cdf_dataset_in.dimensions.items(): cdf_dataset_out.createDimension( name,
    # (len(dimension) if not dimension.isunlimited() else None)) # copy all file data except for the excluded for name,
    # variable in cdf_dataset_in.variables.items(): x = cdf_dataset_out.createVariable(name, variable.datatype,
    # variable.dimensions) cdf_dataset_out[name][:] = cdf_dataset_in[name][:] # copy variable attributes all at once via
    # dictionary cdf_dataset_out[name].setncatts(cdf_dataset_in[name].__dict__)

    # Close the dataset and write the file
    cdf_dataset_out.close()
