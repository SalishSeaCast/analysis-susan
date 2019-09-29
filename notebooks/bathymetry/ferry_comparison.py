import argparse
import netCDF4 as nc
import numpy as np
import pandas as pd
from salishsea_tools import geo_tools, nc_tools, tidetools
import xarray as xr
import glob
import datetime


def get_model_salinity(filenames):
    threemonthsbase = sorted(glob.glob(filenames))
    with nc.Dataset(threemonthsbase[0]) as h:
        model_units = h.variables['time_counter'].units

    with nc_tools.scDataset(threemonthsbase) as f:
        threemonthsbase_sal = f.variables['vosaline'][:, 0, ...]
        timesbase = f.variables['time_counter'][:]
    converted_timesbase = nc.num2date(timesbase, model_units)
    return threemonthsbase_sal, converted_timesbase


def get_pairs(ferry_times, ferry_lons, ferry_lats, ferry_sals,
              ferry_cross, modelfile, ferryfile, accfile, jmodel, imodel, dlon,
              dlat, slon, slat, threemonthsbase_sal, converted_timesbase):
    list_of_modelbase_sals = np.array([])
    list_of_ferrybase_sals = np.array([])
    list_of_lats = np.array([])
    list_of_lons = np.array([])
    list_of_times = np.array([])
    list_of_crossing = np.array([])
    print (threemonthsbase_sal.shape)
    for n in range(ferry_times.shape[0]):
        date = ferry_times[n]
        print (n, date)
        if not ferry_lats.mask[n] and not ferry_sals.mask[n]:
            Xind, Yind = find_point(ferry_lons[n], ferry_lats[n],
                                    jmodel, imodel, dlon, dlat, slon, slat)

            if date.minute <= 30:
                before = (datetime.datetime(year=date.year, month=date.month,
                          day=date.day, hour=date.hour, minute=30) -
                          datetime.timedelta(hours=1))
                index = np.argmin(np.abs(converted_timesbase - date))
                delta = (date - before).seconds / 3600
                s_val = ((delta * (threemonthsbase_sal[index-1, Yind, Xind])) +
                         (1 - delta)*(threemonthsbase_sal[index, Yind, Xind]))
            elif date.minute > 30:
                before = datetime.datetime(year=date.year, month=date.month,
                                           day=date.day, hour=date.hour,
                                           minute=30)
                index = np.argmin(np.abs(converted_timesbase - date))
                delta = (date - before).seconds / 3600
                s_val = ((delta * (threemonthsbase_sal[index, Yind, Xind])) +
                         (1 - delta)*(threemonthsbase_sal[index+1,
                                                          Yind, Xind]))
            list_of_ferrybase_sals = np.append(
                list_of_ferrybase_sals, ferry_sals[n])
            list_of_modelbase_sals = np.append(list_of_modelbase_sals, s_val)
            list_of_lats = np.append(list_of_lats, ferry_lats[n])
            list_of_lons = np.append(list_of_lons, ferry_lons[n])
            list_of_times = np.append(list_of_times, date)
            list_of_crossing = np.append(list_of_crossing, ferry_cross[n])
        if n % 1000 == 0:
            np.savetxt(modelfile, list_of_modelbase_sals)
            np.savetxt(ferryfile, list_of_ferrybase_sals)
    np.savetxt(accfile, list_of_crossing)
    dataout = pd.DataFrame({'time': list_of_times,
                            'lats': list_of_lats,
                            'lons': list_of_lons,
                            'ferry': list_of_ferrybase_sals,
                            'model': list_of_modelbase_sals})
    dataout.to_csv(f'pandas{modelfile}')
    dataout = pd.DataFrame({'crossing': list_of_crossing}, dtype='float')
    dataout.to_csv(f'crossing{modelfile}')
    print(ferry_times[n])
    return


def get_ferry():
    ferry_data = 'ubcONCTWDP1mV18-01_07mar17-02oct17.nc'
#    import ipdb; ipdb.set_trace()
    ferry = nc.Dataset(ferry_data)
    print (ferry)
    unit = ferry.variables['time'].units
    ferry_times = nc.num2date(ferry.variables['time'][:], unit)
    ferry_lats = ferry.variables['latitude'][:]
    ferry_lons = ferry.variables['longitude'][:]
    ferry_sals = ferry.variables['salinity'][:]
    ferry_cross = ferry.variables['crossing_number'][:]
    print ("Start and End Times", ferry_times[0], ferry_times[-1])
    print (ferry_times.shape[0])
    return ferry_times, ferry_lons, ferry_lats, ferry_sals, ferry_cross


def init_find_point():
    with nc.Dataset(
            '/home/sallen/MEOPAR/sea_initial/bathymetry_201803b.nc') as ds:
        nav_lat = ds.variables['nav_lat'][:]
        nav_lon = ds.variables['nav_lon'][:]
    jmodel = np.loadtxt('jvalues.txt', dtype=int)
    imodel = np.loadtxt('ivalues.txt', dtype=int)
    dlon = (nav_lon[imodel[0, -1], jmodel[0, -1]]
            - nav_lon[imodel[0, 0], jmodel[0, 0]]) / jmodel.shape[1]
    dlat = (nav_lat[imodel[-1, 0], jmodel[-1, 0]]
            - nav_lat[imodel[0, 0], jmodel[0, 0]]) / jmodel.shape[0]
    slat = nav_lat[imodel[0, 0], jmodel[0, 0]]
    slon = nav_lon[imodel[0, 0], jmodel[0, 0]]
    return jmodel, imodel, dlon, dlat, slon, slat


def find_point(alon, alat, jmodel, imodel, dlon, dlat, slon, slat):
    i = int((alat - slat + 0.5*dlat)/dlat)
    j = int((alon - slon + 0.5*dlon)/dlon)
    if i >= jmodel.shape[0] or i < 0 or j >= jmodel.shape[1] or j < 0:
        print("find point beyond array", i, j, alat, alon, jmodel.shape)
    return jmodel[i, j], imodel[i, j]


def main(args):
    jmodel, imodel, dlon, dlat, slon, slat = init_find_point()
    print ('done init find point')
    ferry_times,ferry_lons, ferry_lats, ferry_sals, ferry_cross  = get_ferry()
    threemonthsal, converted_timesbase = get_model_salinity(args.datafile)
    get_pairs(ferry_times, ferry_lons, ferry_lats, ferry_sals, ferry_cross,
              args.modelfile, args.ferryfile, args.accfile, jmodel, imodel, dlon, dlat, slon, slat,
              threemonthsal, converted_timesbase)
    print ('Done')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('datafile', help='date files as string with * like: /data/sallen/results/MEOPAR/new_waves/part*/*1h*grid_T*')
    parser.add_argument('ferryfile', help='ferry file name')
    parser.add_argument('modelfile', help='model file name')
    parser.add_argument('accfile', help='extra info file name')
    args = parser.parse_args()
    main(args)
