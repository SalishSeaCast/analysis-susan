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
from pathlib import Path
import sys
import yaml

from salishsea_tools import evaltools as et


def main(config_file):

    with Path(config_file).open("rt") as f:
        config = yaml.safe_load(f)

    start_date = dt.datetime(config['start_date'][0],
                             config['start_date'][1],
                             config['start_date'][2])

    end_date = dt.datetime(config['end_date'][0],
                           config['end_date'][1],
                           config['end_date'][2])

    df1 = et.loadDFOCTD(datelims=(start_date, end_date))
    data = et.matchData(data=df1, filemap=config['filemap'], fdict=config['fdict'],
                        mod_start=start_date, mod_end=end_date,
                        mod_nam_fmt=config['namfmt'], mod_basedir=config['PATH'],
                        mod_flen=config['flen'])
    data.to_csv(config['filename'])


if __name__ == "__main__":
    config_file = sys.argv[1]
    main(config_file)
