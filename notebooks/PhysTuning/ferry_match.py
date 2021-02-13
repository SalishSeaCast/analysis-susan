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
from sqlalchemy import create_engine
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import create_session
from sqlalchemy.sql import and_
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

    dbname = config['sqlfile']
    engine = create_engine('sqlite:///' + dbname)
    Base = automap_base()
    # reflect the tables in salish.sqlite:
    Base.prepare(engine, reflect=True)
    # mapped classes have been created

    # existing tables:
    TDPTurbTBL = Base.classes.TDPTurbTBL
    TDPSalTBL = Base.classes.TDPSalTBL
    GriddedTDPTBL = Base.classes.GriddedTDPTBL

    session = create_session(bind = engine, autocommit = False, autoflush = True)

    df = pd.DataFrame(session.query(GriddedTDPTBL.Year, GriddedTDPTBL.Month,
                                    GriddedTDPTBL.Day, GriddedTDPTBL.Hour,
                                    GriddedTDPTBL.MLat, GriddedTDPTBL.MLon,
                                    GriddedTDPTBL.Modeli.label('i'),
                                    GriddedTDPTBL.Modelj.label('j'),
                                    GriddedTDPTBL.Avg_CT, GriddedTDPTBL.AvgSA,
                                    GriddedTDPTBL.AvgChl_ugl,
                                    GriddedTDPTBL.AvgTurb_NTU
                                    ).filter(and_(GriddedTDPTBL.Year==start_date.year,
                                             GriddedTDPTBL.Month==start_date.month,
                                             GriddedTDPTBL.Modeli>220,
                                             GriddedTDPTBL.Modeli<300,
                                             GriddedTDPTBL.N_SA<10,
                                             GriddedTDPTBL.N_SA>=3)).all())

    df['dtUTC'] = [dt.datetime(iyr,imn,idy,ihr)
                   for ii,(iyr,imn,idy,ihr) in
                   df[['Year','Month','Day','Hour']].iterrows()]

    data=et.matchData(df, filemap=config['filemap'], fdict=config['fdict'],
                      mod_start=start_date, mod_end=end_date,
                      mod_nam_fmt=config['namfmt'], mod_basedir=config['PATH'],
                      mod_flen=config['flen'], method='ferry', preIndexed=True)

    data.to_csv(config['filename'])

    session.close()
    engine.dispose()

if __name__ == "__main__":
    config_file = sys.argv[1]
    main(config_file)
