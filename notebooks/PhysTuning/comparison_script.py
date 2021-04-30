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

"""A collection of Python functions to do routine tasks associated with
plotting and visualization.
"""

import datetime as dt
import pandas as pd
from pathlib import Path
import sys
import yaml

from salishsea_tools import evaltools as et

monthday = {'all': [1, 1, 12, 31],
            '1': [1, 1, 3, 31],
            '2': [3, 1, 6, 30],
            '3': [6, 1, 9, 30],
            '4': [9, 1, 12, 31]}


def main(config_file, obs_data_code, year, quarter):

    start_date = dt.datetime(year, monthday[quarter][0],
                             monthday[quarter][1])
    end_date = dt.datetime(year, monthday[quarter][2],
                           monthday[quarter][3])

    print (start_date, end_date)

    with Path(config_file).open("rt") as f:
        config = yaml.safe_load(f)

    print (obs_data_code, year, quarter)
    preindexed = False
    if obs_data_code == 'ctd':
        df1 = et.loadDFOCTD(config['ctddir'], datelims=(start_date, end_date))
    elif obs_data_code == 'bot':
        df1 = et.loadDFO(basedir=config['botdir'], datelims=(start_date, end_date),
                       excludeSaanich=True)
    elif obs_data_code == 'psf':
        df1 = pd.read_csv(f'{config["csvdir"]}/PSFBotChl.csv', parse_dates=['dtUTC'])
    elif obs_data_code == 'psfts':
        df1 = pd.read_csv(f'{config["csvdir"]}/PSFCTD.csv', parse_dates=['dtUTC'], low_memory=False)
        preindexed = True
    elif obs_data_code == 'pug':
        df1 = pd.read_csv(f'{config["csvdir"]}/WADENuts.csv', parse_dates=['dtUTC'])
    elif obs_data_code == 'pugts':
        df1 = pd.read_csv(f'{config["csvdir"]}/WADECTD.csv', parse_dates=['dtUTC'])
    elif obs_data_code == 'hplc':
        df1 = pd.read_csv(f'{config["csvdir"]}/HPLCPhyto.csv', parse_dates=['dtUTC'])
        preindexed = True
    elif obs_data_code == 'ferry' or obs_data_code == 'ferry_file_only':
        df1 = et.load_ferry_ERDDAP(datelims=(start_date, end_date))
    elif obs_data_code == 'ferry_from_file':
        df1 = pd.read_csv(f'{config["outdir"]}/ferry_{start_date:%Y}.csv', parse_dates=['dtUTC'])
        obs_data_code = 'ferry'
    elif obs_data_code == 'onc' or obs_data_code == 'onc_file_only':
        df1 = et.load_ONC_node_ERDDAP(datelims=(start_date, end_date))
    elif obs_data_code == 'onc_from_file':
        df1 = pd.read_csv(f'{config["outdir"]}/onc_{start_date:%Y}.csv', parse_dates=['dtUTC'])
    else:
        print ('ERROR, specify ctd, bot,  psf, psfts, pug, pugts, hplc, ferry, ferry_from_file, ferry_file_only, onc, onc_from_file, onc_file_only as second argument')

    if obs_data_code == 'ferry_file_only':
        df1.to_csv(f'{config["outdir"]}/ferry_{start_date:%Y}.csv')
    elif obs_data_code == 'onc_file_only':
        df1.to_csv(f'{config["outdir"]}/onc_{start_date:%Y}.csv')
    else:
        if 'ferry' in obs_data_code:
            preindexed = True
            method = 'ferry'
            filename = f'{config["outdir"]}/ObsModel_{config["filebase"]}_ferry_{start_date:%Y%m%d}_{end_date:%Y%m%d}.csv'
        else:
            method = 'bin'
            if 'onc' in obs_data_code:
                preindexed = True
                use_name = 'onc'
            else:
                use_name = obs_data_code
            filename = f'{config["outdir"]}/ObsModel_{config["filebase"]}_{use_name}_{start_date:%Y%m%d}_{end_date:%Y%m%d}.csv'
        data = et.matchData(data=df1, filemap=config['filemap'], fdict=config['fdict'],
                        mod_start=start_date, mod_end=end_date,
                        mod_nam_fmt=config['namfmt'], mod_basedir=config['PATH'],
                        mod_flen=config['flen'], meshPath=config['meshPath'],
                        preIndexed=preindexed, fastSearch=True, method=method)

        print(filename)
        data.to_csv(filename)


if __name__ == "__main__":
    config_file = sys.argv[1]
    obs_data_code = sys.argv[2]
    print(sys.argv)
    year = sys.argv[3]
    quarter = sys.argv[4]
    main(config_file, obs_data_code, int(year), quarter)
