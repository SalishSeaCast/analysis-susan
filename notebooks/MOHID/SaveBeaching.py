# Copyright 2013-- The Salish Sea NEMO Project and
# The University of British Columbia

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# A module to extract beaching times and volumes and
# save boundary with information about the parameters of the spills

import datetime as dt
from pathlib import Path
import xarray as xr

#from StackOverFlow The Aelfinn
def lines_that_contain(string, fp):
    return [line for line in fp if string in line]

def get_parameters(direct):
    for myi in direct.glob('Lagrangian*.dat'):
        asstr = str(myi)
        OilType = asstr[asstr.find('gian_') + 5:asstr.find('-')]
        with open(myi, "r") as fp:
            for line in lines_that_contain("POINT_VOLUME", fp):
                if ':' in line:
                    SpillVolume = float(line[line.find(':')+2:-1])
        with open(myi, "r") as fp:
            for line in lines_that_contain("POSITION_COORDINATES", fp):
                if ':' in line:
                    numbers = (line[line.find(':')+2:-1]).split(' ')
                    lon, lat = float(numbers[0]), float(numbers[1])
    for myi in direct.glob('Model*.dat'):
        with open(myi, "r") as fp:
            for line in lines_that_contain("START", fp):
                numbers = (line[line.find(':')+2:-1]).split(' ')
                startdatetime = dt.datetime(int(numbers[0]), int(numbers[1]), int(numbers[2]),
                                   int(numbers[3]), int(numbers[4]), int(numbers[5]))
    return OilType, SpillVolume, lon, lat, startdatetime


def get_beaching_data(direct):
    for myi in direct.glob('Lagrangian*.nc'):
        data = xr.open_dataset(myi)
        BeachTime = data.Beaching_Time
        BeachVolume = data.Beaching_Volume
        grid_y, grid_x = data.grid_y, data.grid_x

        ncfile = str(myi)
        filename = f'beaching_files/Beaching{(ncfile[ncfile.find("Lagrangian")+10:])}'
    return BeachTime, BeachVolume, grid_y, grid_x, filename


def prepare_dataset(variables, grid_y, grid_x):

    ds_attrs = {
        'acknowledgements':
            'MOHID output',
        'creator_email':
            'sallen@eoas.ubc.ca',
        'creator_name':
            'Salish Sea MEOPAR Project Contributors',
        'creator_url':
            'https://ubc-moad-docs.readthedocs.org/',
        'institution':
            'UBC EOAS',
        'institution_fullname': (
            'Earth, Ocean & Atmospheric Sciences,'
            ' University of British Columbia'
        ),
        'summary': (
            'Beaching Time and Volume from a Specific Run'
        ),
        'source': (
            'analysis-susan/notebooks/MOHID/SaveBeaching.py'
        ),
        'history': (
            '[{}] File creation.'
            .format(dt.datetime.today().strftime('%Y-%m-%d'))
        )
    }

    da = {}
    for var in ['Beaching_Volume', 'Beaching_Time']:
        da[var] = xr.DataArray(
            data=variables[var],
            name=var,
            dims=('grid_y', 'grid_x'),
            coords={
                'grid_y': grid_y,
                'grid_x': grid_x,
            })

    da_attrs = {'OilType': {'units': 'None',
                            'long_name': 'Type of oil spilled and run',
                           },
                'SpillVolume': {'units': 'm3',
                                'long_name': 'Volume of oil initially spilled'},
                'lon': {},
                'lat': {},
                'startdatetime': {'long_name': 'Date and time of Oil Spill'}}
    for var in ['OilType', 'SpillVolume', 'lon', 'lat', 'startdatetime']:
        da[var] = xr.DataArray(
            data=variables[var],
            name=var,
            dims=(),
            attrs=da_attrs[var]
            )

    ds = xr.Dataset(
        data_vars={
            'Beaching_Volume': da['Beaching_Volume'],
            'Beaching_Time': da['Beaching_Time'],
            'OilType': da['OilType'],
            'SpillVolume': da['SpillVolume'],
            'SpillLon': da['lon'],
            'SpillLat': da['lat'],
            'Spilldatetime': da['startdatetime']
        },
        coords={
        'grid_y': grid_y,
                'grid_x': grid_x,
        },
        attrs=ds_attrs
    )

    return ds


def write_out_file(ds, filename):
    encoding = {var: {'zlib': True} for var in ds.data_vars}
    ds.to_netcdf(
        path=filename,
        encoding=encoding,
    )


def SaveBeaching(directory):
    toppath = Path(directory)
    for direct in toppath.glob('*'):
        da = {}
        da['OilType'], da['SpillVolume'], da['lon'], da['lat'], da['startdatetime'] = get_parameters(direct)
        da['Beaching_Time'], da['Beaching_Volume'], grid_y, grid_x, filename = get_beaching_data(direct)
        ds = prepare_dataset(da, grid_y, grid_x)
        write_out_file(ds, filename)


if __name__ == "__main__":
    directory = sys.argv[1]
    SaveBeaching(directory)
