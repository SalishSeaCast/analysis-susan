# Copyright 2013-2016 The Salish Sea MEOPAR contributors
# and The University of British Columbia

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""A Python function to find i,j NEMO coordinates for a fine grid (about 1/10th of a grid point) lat, lon grid
"""

import numpy as np
import sys
import xarray as xr

from salishsea_tools import geo_tools

def main(jstart, jend):

    # load the files
    igrid = np.load(f'igrid_0_{jstart}.npy')
    jgrid = np.load(f'jgrid_0_{jstart}.npy')

    mask_domain = np.load('../mask_domain.npy')

    # calculate the latlon
    decimate_y = 3200
    decimate_x = 1890
    gridcoords = 'coordinates_seagrid_SalishSea201702.nc'
    coords_file = '../../../../grid/'+gridcoords
    coords = xr.open_dataset(coords_file, decode_times=False)

    lon_max = np.floor(coords.glamt.values.max()*decimate_x)/decimate_x
    lon_min = np.floor(coords.glamt.values.min()*decimate_x)/decimate_x
    lat_max = np.floor(coords.gphit.values.max()*decimate_y)/decimate_y
    lat_min = np.floor(coords.gphit.values.min()*decimate_y)/decimate_y
    lons = np.arange(lon_min, lon_max, 1/decimate_x)
    lats = np.arange(lat_min, lat_max, 1/decimate_y)

    # make the loop
    for js, lon in enumerate(lons[jstart:jend]):
        jj = js + jstart
        if int(jj/100)*100 == jj:
            print (jj, lon)
        for ii, lat in enumerate(lats):
            if mask_domain[ii, jj]:
                igrid[ii, jj], jgrid[ii, jj] = geo_tools.find_closest_model_point(
                    lon, lat, coords.glamt[0], coords.gphit[0])

    # write files
    with open(f'igrid_0_{jend}.npy', 'wb') as myfile:
        np.save(myfile, igrid)
    with open(f'jgrid_0_{jend}.npy', 'wb') as myfile:
        np.save(myfile, jgrid)

if __name__ == "__main__":
    jstart = sys.argv[1]
    jend = sys.argv[2]
    main(int(jstart), int(jend))
