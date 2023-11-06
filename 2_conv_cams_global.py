# conda activate CAMS_Globale

# Libreries
import numpy as np
import xarray as xr
import cfgrib
import ecmwflibs
from datetime import datetime, timedelta
import copy
import argparse

from zipfile import ZipFile
import os

def valid_datetime(dt):
    for fmt in ('%Y-%m-%dT%H:%M', '%Y-%m-%dT%H:%M:%S',
                '%Y-%m-%d %H:%M', '%Y-%m-%d %H:%M:%S'):
        try:
            return datetime.strptime(dt, fmt)
        except ValueError:
            pass
    raise argparse.ArgumentTypeError("Invalid date: '{0}'.".format(dt))

def valid_date(d):
    t = 'T00:00'
    return valid_datetime(d + t)

parser = argparse.ArgumentParser(description='Script for converting CAMS Global data set')
parser.add_argument('-s_date', '--start_date', metavar='YYYY-MM-DD HH:MM:SS', type=valid_datetime, required=True)
parser.add_argument('-e_date', '--end_date', metavar='YYYY-MM-DD HH:MM:SS', type=valid_datetime, required=True)
parser.add_argument('-model_order', '--model_order', type=int, required=True)
args = vars(parser.parse_args())

# Temporal interval to consider
# IMPORTANT: the data available of CAMS Global start from 2003-01-02 

# From 2003-1-2 to 2023-1-1
start_datetime_complete = args["start_date"]
end_datetime_complete = args["end_date"]

def joinpath(rootdir, targetdir):
    return os.path.join(os.sep, rootdir + os.sep, targetdir)

# Path of CAMS_Gloabl
path_main_dir_CAMS_Global_data = os.environ['CAMS_Global']

if path_main_dir_CAMS_Global_data == "":
    print("Error: set the environmental variables of CAMS_Global")
    exit(-1)

DATADIR = joinpath(path_main_dir_CAMS_Global_data, "datasets_model_level_" + str(model_order))

DATADIR_AIR_POLL_ORG = joinpath(DATADIR, "EAC4_ORG")
DATADIR_RELATIVE_HUMIDITY= joinpath(DATADIR, "relative_humidity")

DATADIR_ALL_AIR_POLL_FORMULA = joinpath(DATADIR, "italy_ext_formula")

if not os.path.exists(DATADIR_ALL_AIR_POLL_FORMULA):
    os.mkdir(DATADIR_ALL_AIR_POLL_FORMULA)

DATADIR_ALL_AIR_POLL_FIXED = joinpath(DATADIR, "italy_ext_fixed")

if not os.path.exists(DATADIR_ALL_AIR_POLL_FIXED):
    os.mkdir(DATADIR_ALL_AIR_POLL_FIXED)

# IMPORTANT: the ozone measure O3 is saved as GO3

# List of all air pollutants (Mass Mixing ratio)
# NO2, SO2, CO, O3: kg kg^-1
# PM1, PM2.5, PM10: kg m^-3
list_air_poll = ["NO2", "GO3", "SO2", "CO"]
list_air_poll_pm = ["PM1", "PM2p5", "PM10"]

# ["NO2", "GO3", "SO2", "CO"]
for air_chem in list_air_poll:
    
    DATADIR_CURRENT_AIR_POLL_FORMULA = joinpath(DATADIR_ALL_AIR_POLL_FORMULA, air_chem)
    DATADIR_CURRENT_AIR_POLL_FIXED = joinpath(DATADIR_ALL_AIR_POLL_FIXED, air_chem)

    if not os.path.exists(DATADIR_CURRENT_AIR_POLL_FORMULA):
        os.mkdir(DATADIR_CURRENT_AIR_POLL_FORMULA)
    
    if not os.path.exists(DATADIR_CURRENT_AIR_POLL_FIXED):
        os.mkdir(DATADIR_CURRENT_AIR_POLL_FIXED)

# ["PM1", "PM2p5", "PM10"]
for air_chem in list_air_poll_pm:
    
    DATADIR_CURRENT_AIR_POLL_FORMULA = joinpath(DATADIR_ALL_AIR_POLL_FORMULA, air_chem)
    DATADIR_CURRENT_AIR_POLL_FIXED = joinpath(DATADIR_ALL_AIR_POLL_FIXED, air_chem)

    if not os.path.exists(DATADIR_CURRENT_AIR_POLL_FORMULA):
        os.mkdir(DATADIR_CURRENT_AIR_POLL_FORMULA)
    
    if not os.path.exists(DATADIR_CURRENT_AIR_POLL_FIXED):
        os.mkdir(DATADIR_CURRENT_AIR_POLL_FIXED)

def E_s_T(temp):
  
    eso = 6.1078

    c4 = 0.43884187 * 1e-8

    c = [   
            0.99999683, -0.90826951*1e-2, 0.78736169*1e-4, -0.61117958*1e-6, \
            c4, -0.29883885*1e-10, 0.21874425*1e-12, -0.17892321*1e-14, \
            0.11112018*1e-16, -0.30994571*1e-19
        ]

    p = (c[0]+temp*(c[1]+temp*(c[2]+temp* \
        (c[3]+temp*(c[4]+temp*(c[5]+temp*(c[6]+temp* \
        (c[7]+temp*(c[8]+temp*(c[9]))))))))))

    # Polynomial function
    E_s = eso / pow(p,8)

    return E_s

def compute_air_density(ds_chem, ds_chem_rh):

    # DEFINITIONS OF ATMOSPHERE LEVELS OF CAMS:
    # https://confluence.ecmwf.int/display/UDOC/L60+model+level+definitions

    # Temperature in Kelvin
    temperature_ds = ds_chem['t']

    # Temperature in Celsius
    temperature_ds_C = temperature_ds - 273.15

    # Relative umidity
    rh_ds = ds_chem_rh['r']

    # Surface pressure expressed in Pa 
    # (atmosphere pressure of model level 55 --> about 288 metri = 985,15 hPa)
    surface_pressure_Pa = 1013.25 * 100

    # Specific constant of dry air: 287.05 J/(kgK)
    constant_dry_air = 287.05

    # Specific constant of water vapor: 461.495 J/(kgK)
    constant_h2o = 461.495

    P_v = rh_ds * E_s_T(temperature_ds_C)
    
    # Compute the density of water vapor (kg/m^3)
    h2o_density = P_v / (constant_h2o * temperature_ds)

    P_d = surface_pressure_Pa - P_v

    # Compute the density of dry air (kg/m^3)
    # Formula: Atmosphere_pressure(hPa) /(constant_dry_air(J/(khK)) * Air_temp(K))
    dry_air_density = P_d / (constant_dry_air * temperature_ds)

    # Compute air density (kg/m^3)
    air_density = dry_air_density + h2o_density

    return air_density

def conv_to_concentration(ds_chem, current_start_date, delta, list_of_air_pollutants, diff_dates_hours, air_density, is_fixed):

    for air_chem in list_of_air_pollutants:
        
        current_date_loc = current_start_date

        for time in range(diff_dates_hours):
            
            current_ds_chem = copy.copy(ds_chem[air_chem.lower()].sel(time=current_date_loc.isoformat()))

            if air_chem == "NO2" or air_chem == "SO2" or air_chem == "GO3" or air_chem == "CO":
        
                # -------- Conversion from MMR (kg kg^-1) --> concentration (kg m^-3) --------
                if is_fixed:
                    current_air_poll_ds = current_ds_chem * air_density
                else:
                    current_air_poll_ds = current_ds_chem * air_density[time,:,:]

                if air_chem == "NO2" or air_chem == "SO2" or air_chem == "GO3":
    
                    # -------- Conversion (kg m^-3) --> (ug m^-3) --------
                    current_air_poll_ds *= 1e+9

                else:
            
                    # -------- Conversion (kg m^-3) --> (mg m^-3) --------
                    current_air_poll_ds *=  1e+6
            
                ds_chem[air_chem.lower()].loc[dict(time=current_date_loc.isoformat())] = current_air_poll_ds

            elif air_chem == "PM1" or air_chem == "PM2p5" or air_chem == "PM10":
        
                # -------- Conversion (kg m^-3) --> (ug m^-3) --------
                current_air_poll_ds = current_ds_chem * 1e+9

                ds_chem[air_chem.lower()].loc[dict(time=current_date_loc.isoformat())] = current_air_poll_ds

            current_date_loc += delta

    # Remove the last observation of the first day of next month.
    ds_chem = ds_chem.sel(time=slice(current_start_date, current_date_loc - delta))

    return ds_chem

def func_ext_save_air_poll(ds_chem, list_to_consider, main_dir, string_current_date, is_fixed):

    for air_chem in list_to_consider:
            
        DATADIR_CURRENT_AIR_POLL = joinpath(main_dir, air_chem)

        DATADIR_CURRENT_AIR_POLL = joinpath(DATADIR_CURRENT_AIR_POLL, string_current_date)

        if not os.path.exists(DATADIR_CURRENT_AIR_POLL):
            os.mkdir(DATADIR_CURRENT_AIR_POLL)
            
        path_netcdf_file = joinpath(DATADIR_CURRENT_AIR_POLL, string_current_date + ".nc")

        ds_air_chem = ds_chem[air_chem.lower()]
        ds_air_chem.to_netcdf(path_netcdf_file)

        if is_fixed:
            print(air_chem + " FIXED " + string_current_date + " analysed")
        else:
            print(air_chem + " FORMULA " + string_current_date + " analysed")

# Time resolution of CAMS Global
time_res = 3

delta = timedelta(hours=time_res)

esito = False
current_date = start_datetime_complete

while esito == False:

    # Andiamo a scalare di un giorno al fine di definire l'ultimo giorno del mese
    if current_date.date() > end_datetime_complete.date():
        current_date = datetime(end_datetime_complete.year, end_datetime_complete.month, 1, 0, 0)
        end_current_date = end_datetime_complete
        esito = True
    else:
        if current_date.month + 1 == 13:
            end_current_date = datetime(current_date.year + 1, 1, 1, 0, 0)
        else:
            end_current_date = datetime(current_date.year, current_date.month + 1, 1, 0, 0)
    
    # Loading data sets of ["NO2", "GO3", "SO2", "CO"]
    name_file_air_poll_netcdf = "air_poll_" + str(current_date.year) + "-" + str(current_date.month).zfill(2) + ".grib"
    PATH_NETCDF_FILE_AIR_POLL = joinpath(DATADIR_AIR_POLL_ORG, name_file_air_poll_netcdf)

    # Loading data sets of ["PM1", "PM2p5", "PM10"]
    name_file_pm_air_poll_netcdf = "PM_air_poll_" + str(current_date.year) + "-" + str(current_date.month).zfill(2) + ".grib"
    PATH_NETCDF_FILE_PM_AIR_POLL = joinpath(DATADIR_AIR_POLL_ORG, name_file_pm_air_poll_netcdf)

    # Loading Relative umidity
    name_file_relative_humidity_netcdf = "relative_humidity_" + str(current_date.year) + "-" + str(current_date.month).zfill(2) + ".nc"
    PATH_NETCDF_FILE_RELATIVE_HUMIDITY = joinpath(DATADIR_RELATIVE_HUMIDITY, name_file_relative_humidity_netcdf)

    ds_chem_formula = cfgrib.open_dataset(PATH_NETCDF_FILE_AIR_POLL)
    ds_chem_fixed = cfgrib.open_dataset(PATH_NETCDF_FILE_AIR_POLL)

    ds_chem_formula_pm = cfgrib.open_dataset(PATH_NETCDF_FILE_PM_AIR_POLL)
    ds_chem_fixed_pm = cfgrib.open_dataset(PATH_NETCDF_FILE_PM_AIR_POLL)

    ds_chem_rh = xr.open_dataset(PATH_NETCDF_FILE_RELATIVE_HUMIDITY)

    idx_lat = ds_chem_formula.indexes["latitude"]
    idx_lon = ds_chem_formula.indexes["longitude"]

    # Compute the air density
    air_density_formula = compute_air_density(ds_chem_formula, ds_chem_rh)

    # Definition of air density fixed of model level 55 --> about 288 meters
    # https://confluence.ecmwf.int/display/UDOC/L60+model+level+definitions
    air_density_fixed = 1.191403

    diff_dates = end_current_date - current_date
    diff_dates_hours = int(diff_dates.total_seconds() / ((60*60)*time_res))

    # -------------------- Conversion da MMR (kg kg^-1) --> concentration (kg m^-3) --------------------
    # Formula: MMR * Air density --> (kg kg^-1) * (kg m^-3) --> (kg m^-3)
    ds_chem_formula = conv_to_concentration(ds_chem_formula, current_date, delta, list_air_poll, diff_dates_hours, air_density_formula, False)
    ds_chem_fixed = conv_to_concentration(ds_chem_fixed, current_date, delta, list_air_poll, diff_dates_hours, air_density_fixed, True)

    ds_chem_formula_pm = conv_to_concentration(ds_chem_formula_pm, current_date, delta, list_air_poll_pm, diff_dates_hours, air_density_formula, False)
    ds_chem_fixed_pm = conv_to_concentration(ds_chem_fixed_pm, current_date, delta, list_air_poll_pm, diff_dates_hours, air_density_fixed, True)

    string_current_date = str(current_date.year) + "-" + str(current_date.month).zfill(2)
    
    func_ext_save_air_poll(ds_chem_formula, list_air_poll, DATADIR_ALL_AIR_POLL_FORMULA, string_current_date, False)
    func_ext_save_air_poll(ds_chem_fixed, list_air_poll, DATADIR_ALL_AIR_POLL_FIXED, string_current_date, True)

    func_ext_save_air_poll(ds_chem_formula_pm, list_air_poll_pm, DATADIR_ALL_AIR_POLL_FORMULA, string_current_date, False)
    func_ext_save_air_poll(ds_chem_fixed_pm, list_air_poll_pm, DATADIR_ALL_AIR_POLL_FIXED, string_current_date, True)

    current_date = end_current_date