# conda activate CAMS_Global

# Libreries
import numpy as np
import xarray as xr
import cfgrib
import ecmwflibs
from zipfile import ZipFile
import os
import argparse
import warnings

# API Request
import cdsapi
from datetime import datetime, timedelta

warnings.filterwarnings("ignore")

def joinpath(rootdir, targetdir):
    return os.path.join(os.sep, rootdir + os.sep, targetdir)

# Path of CAMS_Gloabl
path_main_dir_CAMS_Global_data = os.environ['CAMS_Global']

if path_main_dir_CAMS_Global_data == "":
    print("Error: set the environmental variables of CAMS_Global")
    exit(-1)

if not os.path.exists(path_main_dir_CAMS_Global_data):
  os.mkdir(path_main_dir_CAMS_Global_data)

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

parser = argparse.ArgumentParser(description='Script for downloading CAMS Global data set')
parser.add_argument('-s_date', '--start_date', metavar='YYYY-MM-DD HH:MM:SS', type=valid_datetime, required=True)
parser.add_argument('-e_date', '--end_date', metavar='YYYY-MM-DD HH:MM:SS', type=valid_datetime, required=True)
parser.add_argument('-model_level', '--model_level', type=int, required=True)
args = vars(parser.parse_args())

# Temporal interval to consider
start_datetime_complete = args["start_date"]
end_datetime_complete = args["end_date"]

model_level = int(args["model_level"])

# Main directory of CAMS Global where will be saved the NetCDF files 
DATADIR = joinpath(path_main_dir_CAMS_Global_data, "datasets_model_level_" + str(model_level))

if not os.path.exists(DATADIR):
  os.mkdir(DATADIR)

DATADIR_AIR_POLL_ORG = joinpath(DATADIR, "EAC4_ORG")

if not os.path.exists(DATADIR_AIR_POLL_ORG):
  os.mkdir(DATADIR_AIR_POLL_ORG)

DATADIR_RELATIVE_HUMIDITY= joinpath(DATADIR, "relative_humidity")

if not os.path.exists(DATADIR_RELATIVE_HUMIDITY):
  os.mkdir(DATADIR_RELATIVE_HUMIDITY)

# Coordinate dell'Italia
lat_italy_bnds, lon_italy_bnds = [32,50], [5,21]

c = cdsapi.Client(url='https://ads.atmosphere.copernicus.eu/api/v2', key='15833:3671670b-6ec7-47c5-9483-8175ee8648c9')

# List of all air pollutants (Mass Mixing ratio)
# NO2, SO2, CO, O3: kg kg^-1
# PM1, PM2.5, PM10: kg m^-3
list_of_air_pollutants = ["pm1", "pm2p5", "pm10", "co", "no2", "o3", "so2"]

# Time resolution of CAMS Global
time_res = 3

delta = timedelta(hours=time_res)

current_date = start_datetime_complete.date()
esito = False

while esito == False:

    if current_date > end_datetime_complete.date():
        current_date = datetime(end_datetime_complete.year, end_datetime_complete.month, 1, 0, 0).date()
        end_current_date = end_datetime_complete.date()
        esito = True
    else:
        if current_date.month + 1 == 13:
            end_current_date = datetime(current_date.year + 1, 1, 1, 0, 0).date()
        else:
            end_current_date = datetime(current_date.year, current_date.month + 1, 1, 0, 0).date()
    
    name_file_air_poll_netcdf = "air_poll_" + str(current_date.year) + "-" + str(current_date.month).zfill(2) + ".grib"
    PATH_NETCDF_FILE_AIR_POLL = joinpath(DATADIR_AIR_POLL_ORG, name_file_air_poll_netcdf)

    if not os.path.exists(PATH_NETCDF_FILE_AIR_POLL):

        c.retrieve(
            'cams-global-reanalysis-eac4',
            {
                'format': 'grib',
                'variable': [
                    'carbon_monoxide', 'nitrogen_dioxide',
                    'ozone', 'sulphur_dioxide', 'temperature',
                ],
                'model_level': str(model_level),
                'time': [
                    '00:00', '03:00', '06:00',
                    '09:00', '12:00', '15:00',
                    '18:00', '21:00',
                ],
                'area': [lat_italy_bnds[1], lon_italy_bnds[0], lat_italy_bnds[0], lon_italy_bnds[1]],
                'date': current_date.isoformat() + "/" + end_current_date.isoformat(),
            },
            PATH_NETCDF_FILE_AIR_POLL)

        print("Orginal data of EAC4 (NO2-O3-SO2-CO) " + str(current_date.year) + "-" + str(current_date.month).zfill(2) + " converted is saved")

    name_file_pm_air_poll_netcdf = "PM_air_poll_" + str(current_date.year) + "-" + str(current_date.month).zfill(2) + ".grib"
    PATH_NETCDF_FILE_PM_AIR_POLL = joinpath(DATADIR_AIR_POLL_ORG, name_file_pm_air_poll_netcdf)

    if not os.path.exists(PATH_NETCDF_FILE_PM_AIR_POLL):
        c.retrieve(
            'cams-global-reanalysis-eac4',
            {
                'format': 'grib',
                'variable': [
                    'particulate_matter_10um', 'particulate_matter_1um',
                    'particulate_matter_2.5um',
                ],

                'time': [
                    '00:00', '03:00', '06:00',
                    '09:00', '12:00', '15:00',
                    '18:00', '21:00',
                ],
                'area': [ lat_italy_bnds[1], lon_italy_bnds[0], lat_italy_bnds[0], lon_italy_bnds[1]],
                'date': current_date.isoformat() + "/" + end_current_date.isoformat(),
            },
            PATH_NETCDF_FILE_PM_AIR_POLL)

        print("Orginal data of EAC4 (PM_1-PM2,5-PM10) " + str(current_date.year) + "-" + str(current_date.month).zfill(2) + " converted is saved")
    
    name_file_relative_humidity_netcdf = "relative_humidity_" + str(current_date.year) + "-" + str(current_date.month).zfill(2) + ".nc"
    PATH_NETCDF_FILE_RELATIVE_HUMIDITY = joinpath(DATADIR_RELATIVE_HUMIDITY, name_file_relative_humidity_netcdf)

    if not os.path.exists(PATH_NETCDF_FILE_RELATIVE_HUMIDITY):

        c.retrieve(
            'cams-global-reanalysis-eac4',
            {
                'format': 'netcdf',
                'pressure_level': '1000',
                'variable': 'relative_humidity',
                'time': [
                    '00:00', '03:00', '06:00',
                    '09:00', '12:00', '15:00',
                    '18:00', '21:00',
                ],
                'area': [ lat_italy_bnds[1], lon_italy_bnds[0], lat_italy_bnds[0], lon_italy_bnds[1]],
                'date': current_date.isoformat() + "/" + end_current_date.isoformat(),
            },
            PATH_NETCDF_FILE_RELATIVE_HUMIDITY)

    print("Relative humidity " + str(current_date.year) + "-" + str(current_date.month).zfill(2) + " converted is saved")

    current_date = end_current_date