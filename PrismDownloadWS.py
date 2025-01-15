"""
Anar Batmunkh

PRISM Web Service Documentation: https://prism.oregonstate.edu/documents/PRISM_downloads_web_service.pdf

DATA AVAILABILITY

Currently, the following data (i.e., grids in BIL format) can be downloaded via the web service:
    • Variables:
        o ppt, tmin, tmax, tmean, tdmean, vpdmin, vpdmax (as both time series grids and normals grids)
        o solslope, soltotal, solclear, soltrans (as normals grids only)
    • Recent data: 1981 to present; daily, monthly, annual data
    • Historical data: 1895 through 1980; complete monthly and annual data by year
    • 30-year normals data: daily, monthly, and annual normals, each as a single grid

DOWNLOAD LIMITS

Download activity is continuously monitored. To prevent rogue download scripts from exceeding bandwidth limits,
if a file is downloaded twice in a 24-hour period, no more downloads of that file will be allowed during that period.
Repeated excessive download activity may result in IP address blocking, at our discretion.

DATA WEB SERVICE USAGE
Time Series (daily, monthly, annual) Data Request Syntax

To initiate a request via this web service, use a client (e.g., web browser, the wget or curl utility, etc.) that can
perform an HTTP request with the following parameters:

    https://services.nacse.org/prism/data/public/4km/<element>/<date><?format=[nc|asc|grib2]>
        Where:
            • <element> is ppt, tmin, tmax, tmean, tdmean, vpdmin, or vpdmax
            • <date> is:
                o YYYYMMDD for daily data (between yesterday and January 1st, 1981) – returns a single grid in a
                  .zip file
                o YYYYMM for monthly data (between last month and January 1981) – returns a single grid in a .zip
                  file
                o YYYY for annual data (between last year and 1981) – returns a single grid in a .zip file
                o YYYY for historical data (between 1980 and 1895) – returns a single zip file containing 12 monthly
                  grids for YYYY plus the annual.
            • <?format=[nc|asc|grib2]>: optional command that will deliver the data in the specified format.
              Format options:
                o nc: netCDF format (https://www.unidata.ucar.edu/software/netcdf)
                o asc: ASCII Grid format (https://www.loc.gov/preservation/digital/formats/fdd/fdd000421.shtml)
                o grib2: Grib2 format (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc)
                
Valid HTTP examples:

Daily grid: https://services.nacse.org/prism/data/public/4km/tmin/20090405
Monthly grid: https://services.nacse.org/prism/data/public/4km/ppt/198304
Annual grid: https://services.nacse.org/prism/data/public/4km/tmax/2021
Historical monthly grid: https://services.nacse.org/prism/data/public/4km/tmin/194402
Reformatted (netCDF) monthly grid: https://services.nacse.org/prism/data/public/4km/tmin/200904?format=nc
Reformatted (Grib2) daily grid: https://services.nacse.org/prism/data/public/4km/ppt/19911023?format=grib2

"""


import os
import pandas as pd
import numpy as np
import datetime
import urllib.request
import zipfile
from osgeo import gdal
from tqdm import tqdm
import time


gdal.UseExceptions()


def get_value_v2(rasterfile, points):
    gdata = gdal.Open(rasterfile)
    gt = gdata.GetGeoTransform()
    data = gdata.ReadAsArray().astype(np.float64)
    gdata = None
    out_data = []
    for i in range(len(points)):
        pos = (points[i][1], points[i][2])  # Using the order longitude, latitude
        x = int((pos[0] - gt[0]) / gt[1])
        y = int((pos[1] - gt[3]) / gt[5])
        try:
            value = data[y, x]
        except IndexError:
            value = np.nan  # Handle case where the index is out-of-range
        out_data.append(value)
    return np.array(out_data)


def download_prism_data(url, file_name, data_path):
    os.chdir(data_path)
    urllib.request.urlretrieve(url, file_name)
    
    try:
        with zipfile.ZipFile(file_name, 'r') as zip_ref:
            zip_ref.extractall()
    except zipfile.BadZipFile:
        print(f"Downloaded file {file_name} is not a valid ZIP file.")
    finally:
        if os.path.isfile(file_name):
            os.remove(file_name)


def extract_prism_data(data_path, target_day, data_type):
    files = os.listdir(data_path)
    target_files = [f for f in files if f.endswith('.bil') and data_type in f]
    target_files.sort(reverse=True)
    src_filename = None
    for filename in target_files:
        if filename[-16:-8] == target_day:
            src_filename = filename
            break
    return src_filename


# Input CSV with Field ID and coordinates (Field ID, longitute, latitude)
input_CSV = r"U:\z_CornOnFarmBioN\Anar_Tools\2023-2024 mizzou field location.csv"
coordinates_df = pd.read_csv(input_CSV)
points = coordinates_df.to_numpy()

# Date in range
start_date = datetime.date(2017, 5, 1)  # Start year
end_date = datetime.date(2019, 5, 1)  # End year

# Temporary datasets downloaded from the web service will be stored here
data_path = r"U:\z_CornOnFarmBioN\Anar_Tools\PRISM_data"

# CSV files will be stored here (Field ID\CSVs)
output_base_path = r"U:\z_CornOnFarmBioN\Anar_Tools\Output_CSV"

os.makedirs(data_path, exist_ok=True)
os.makedirs(output_base_path, exist_ok=True)

# Start the loop from the YEAR
for year in range(start_date.year, end_date.year + 1):
    if year == start_date.year:
        year_start_date = start_date
    else:
        year_start_date = datetime.date(year, 1, 1)
    
    if year == end_date.year:
        year_end_date = end_date
    else:
        year_end_date = datetime.date(year, 12, 31)

    # Data structure to hold data for the current year
    yearly_data = {field_id: [] for field_id in coordinates_df['Field ID'].unique()}

    for single_date in tqdm(pd.date_range(year_start_date, year_end_date)):
        date_str = single_date.strftime('%Y%m%d')
        try:
            # Pass the URLs for different data types (see PRISM document below)
            # https://prism.oregonstate.edu/documents/PRISM_downloads_web_service.pdf
            base_url = 'http://services.nacse.org/prism/data/public/4km'
            urls = {
                'ppt': f'{base_url}/ppt/{date_str}',
                'tmin': f'{base_url}/tmin/{date_str}',
                'tmax': f'{base_url}/tmax/{date_str}',
                'tmean': f'{base_url}/tmean/{date_str}'
            }
            

            values = {}
            for key, url in urls.items():
                file_name = f'{date_str}_{key}.zip'
                download_prism_data(url, file_name, data_path)

                time.sleep(1)
                
                src_filename = extract_prism_data(data_path, date_str, key)
                if src_filename:
                    values[key] = get_value_v2(src_filename, points)
                else:
                    values[key] = np.full(len(points), np.nan)

            # Collect data for each coordinate (XY)
            for idx, row in coordinates_df.iterrows():
                field_id = row['Field ID']
                yearly_data[field_id].append([
                    row['Field ID'], row['longitude'], row['latitude'], date_str,
                    values['ppt'][idx], values['tmin'][idx], values['tmax'][idx], values['tmean'][idx]
                ])

            # Delete the extracted files from the download
            for file in os.listdir(data_path):
                file_path = os.path.join(data_path, file)
                if os.path.isfile(file_path): 
                    os.remove(file_path)

        except Exception as e:
            print(f"Error processing date {date_str}: {e}")

    # Save the yearly data to CSV files
    for field_id, data in yearly_data.items():
        field_dir = os.path.join(output_base_path, str(field_id))
        if not os.path.exists(field_dir):
            os.makedirs(field_dir)

        output_df = pd.DataFrame(data, columns=['Field ID', 'longitude', 'latitude', 'date', 'precipitation', 'min_temp', 'max_temp', 'mean_temp'])
        output_file = os.path.join(field_dir, f"{year}.csv")
        output_df.to_csv(output_file, index=False)



