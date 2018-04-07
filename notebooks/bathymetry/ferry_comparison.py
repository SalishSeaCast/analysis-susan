import argparse
import netCDF4 as nc
import numpy as np
from salishsea_tools import geo_tools, nc_tools, tidetools
import xarray as xr
import glob
import datetime

def get_model_salinity(filenames):
    threemonthsbase = sorted(glob.glob(filenames))
    with nc.Dataset(threemonthsbase[0]) as h:
        model_units = h.variables['time_counter'].units

    with nc_tools.scDataset(threemonthsbase) as f:
        threemonthsbase_sal = f.variables['vosaline'][:,1,...]
        timesbase = f.variables['time_counter'][:]
    converted_timesbase = nc.num2date(timesbase, model_units)
    return threemonthsbase_sal, converted_timesbase


def get_pairs(istart, iend, ferry_times, ferry, modelfile, ferryfile, accfile, jmodel, imodel, dlon, dlat, slon, slat,
              threemonthsbase_sal, converted_timesbase):
    list_of_modelbase_sals = np.array([])
    list_of_ferrybase_sals = np.array([])
    list_of_lats = np.array([])
    list_of_lons = np.array([])
    list_of_times = np.array([])
    list_of_crossing = np.array([])
    for n in range(istart, iend):
        print ('start')
        date = ferry_times[n]
        if ((ferry.variables['s.latitude'][n].mask == False)
                 and (ferry.variables['s.salinity'][n].mask == False)):
            print ('switch')
            Xind, Yind = find_point(ferry.variables['s.longitude'][n], ferry.variables['s.latitude'][n],
                                    jmodel, imodel, dlon, dlat, slon, slat)
            print ('find')

            if date.minute <= 30:
                before = datetime.datetime(year = date.year, month = date.month, day = date.day,
                                           hour = (date.hour), minute = 30) - datetime.timedelta(hours=1)
                index = np.argmin(np.abs(converted_timesbase - date))
                delta = (date - before).seconds / 3600
                s_val = ((delta * (threemonthsbase_sal[index-1, Yind, Xind])) +
                         (1- delta)*(threemonthsbase_sal[index, Yind, Xind]))
            elif date.minute > 30:
                before = datetime.datetime(year = date.year, month = date.month, day = date.day,
                                           hour = (date.hour), minute = 30)
                index = np.argmin(np.abs(converted_timesbase - date))
                delta = (date - before).seconds / 3600
                s_val = ((delta * (threemonthsbase_sal[index, Yind, Xind])) +
                         (1- delta)*(threemonthsbase_sal[index+1, Yind, Xind]))
            print('calculate')
            list_of_ferrybase_sals = np.append(list_of_ferrybase_sals, ferry.variables['s.salinity'][n])
            list_of_modelbase_sals = np.append(list_of_modelbase_sals, s_val)
            list_of_lats = np.append(list_of_lats, ferry.variables['s.latitude'][n])
            list_of_lons = np.append(list_of_lons, ferry.variables['s.longitude'][n])
            list_of_crossing = np.append(list_of_crossing, ferry.variables['s.crossing_number'][n])
            print(n)
        if n % 1000 == 0:
            np.savetxt(modelfile, list_of_modelbase_sals)
            np.savetxt(ferryfile, list_of_ferrybase_sals)
            np.savetxt(accfile, (list_of_lats, list_of_lons, list_of_crossing))
            print(ferry_times[n])
    np.savetxt(modelfile, list_of_modelbase_sals)
    np.savetxt(ferryfile, list_of_ferrybase_sals)
    np.savetxt(accfile, (list_of_lats, list_of_lons, list_of_crossing))
    print(ferry_times[n])
    return


def get_ferry(istart, iend):
    ferry_data = 'https://salishsea.eos.ubc.ca/erddap/tabledap/ubcONCTWDP1mV18-01'
    ferry = nc.Dataset(ferry_data)
    unit = ferry.variables['s.time'].units
    ferry_times = nc.num2date(ferry.variables['s.time'][:], unit)
    print ("Start and End Times", ferry_times[istart], ferry_times[iend])
    return ferry, ferry_times, unit


def init_find_point():
    with nc.Dataset('/home/sallen/MEOPAR/sea_initial/bathymetry_201803b.nc') as ds:
         nav_lat = ds.variables['nav_lat'][:]
         nav_lon = ds.variables['nav_lon'][:]
    jmodel = np.loadtxt('jvalues.txt', dtype=int)
    imodel = np.loadtxt('ivalues.txt', dtype=int)
    dlon = (nav_lon[imodel[0, -1], jmodel[0, -1]] - nav_lon[imodel[0, 0], jmodel[0, 0]]) / jmodel.shape[1]
    dlat = (nav_lat[imodel[-1, 0], jmodel[-1, 0]] - nav_lat[imodel[0, 0], jmodel[0, 0]]) / jmodel.shape[0]
    slat = nav_lat[imodel[0,0], jmodel[0,0]]
    slon = nav_lon[imodel[0,0], jmodel[0,0]]
    return jmodel, imodel, dlon, dlat, slon, slat


def find_point(alon, alat, jmodel, imodel, dlon, dlat, slon, slat):
    i = int((alat - slat + 0.5*dlat)/dlat)
    j = int((alon - slon + 0.5*dlon)/dlon)
    return jmodel[i, j], imodel[i, j]


def main(args):
    jmodel, imodel, dlon, dlat, slon, slat = init_find_point()
    ferry, ferry_times, unit = get_ferry(int(args.istart), int(args.iend))
    threemonthsal, converted_timesbase = get_model_salinity(args.datafile)
    get_pairs(int(args.istart), int(args.iend), ferry_times, ferry,
              args.modelfile, args.ferryfile, args.accfile, jmodel, imodel, dlon, dlat, slon, slat,
              threemonthsal, converted_timesbase)
    print ('Done')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('datafile', help='date files as string with * like: /data/sallen/results/MEOPAR/new_waves/part*/*1h*grid_T*')
    parser.add_argument('ferryfile', help='ferry file name')
    parser.add_argument('modelfile', help='model file name')
    parser.add_argument('accfile', help='extra info file name')
    parser.add_argument('istart', help='first index')
    parser.add_argument('iend', help='last index')
    args = parser.parse_args()
    main(args)

