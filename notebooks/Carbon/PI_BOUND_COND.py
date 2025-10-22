import datetime
import gsw
import numpy as np
import pandas as pd
import PyCO2SYS as pyco2
import xarray as xr


# Calculate the Atmospheric CO2 when watermass was at surface: atpco2
def calculate_atpco2_for_years(years, include_seasonal=False):
    co2_rec = pd.read_csv('../../../analysis_tereza/notebooks/PI_CARBON_PAPER/PI_BOUND_COND/CLEAN/lawdome_maunaloa.csv')
    zz_year = years.min()
    if (zz_year > 1832) & (zz_year < 1950): # Note, just using one value here!
        atpco2 = co2_rec.PPMCO2[co2_rec.YEAR == int(zz_year)].item()
    elif zz_year < 2000: # Use the parabola fit
        sindex = co2_rec[co2_rec['YEAR'] == 1950].index.values[0]
        co2_rec = co2_rec.set_index('YEAR')
        c = np.polynomial.Polynomial.fit(co2_rec.index[sindex:], co2_rec.PPMCO2[sindex:], deg=2)
        atpco2 = c(years)
    elif zz_year < 2025: # Use SalishSeaCast linear fit
        zz_LR_slope   =  2.149      
        zz_LR_int     =  -3929.359
        atpco2 = years * zz_LR_slope + zz_LR_int 
    else:
        print ('Year asked for is not in range')
        stop
    if include_seasonal:
        zz_ctr        =  161.898   #  based on scripps observations at ptbarrow and lajolla
        zz_amp        =  7.083     #  calc notebook by TJSJ
        zz_wid        =  44.703    #  notebook loc:
        zz_ctr2       =  218.832   #  carbon_dev/MOCSY_and_FLUX/CO2_obs.ipynb
        zz_amp2       =  -19.004   #
        zz_wid2       =  87.8836   #
        zz_ctr3       =  199.430   #
        zz_amp3       =  8.026     #
        zz_wid3       =  -185.920  #
        zz_days = years - np.trunc(years)
        zz_yearcyc = ( zz_amp * np.exp( -((zz_days - zz_ctr)/zz_wid)**2) 
              + zz_amp2 * np. exp( -((zz_days - zz_ctr2)/zz_wid2)**2) 
              + zz_amp3 * np.exp( -((zz_days - zz_ctr3)/zz_wid3)**2)
             )
        atpco2 = atpco2 + zz_yearcyc
    return atpco2
    
# Calculate the Atmospheric CO2 when watermass was at surface: atpco2
def calculate_atpco2_for_watercolumn(ds):
    atpco2 = calculate_atpco2_for_years(ds.pycnal_last_at_surface.values)
    ds['atpco2'] = xr.DataArray(atpco2,
                                  coords=ds.coords,
                               attrs={'long_name': 'Atmospheric pCO2',
                                        'units': 'ppm'})
    return ds


# Calculate surface pCO2 for now and for pre-industrial
def surface_atpco2_for_year(ds, now, preindustrial, include_seasonal=False):
    atpco2 = calculate_atpco2_for_years(np.array([now])), calculate_atpco2_for_years(np.array([preindustrial]), include_seasonal)
    return atpco2
    

# Calculate Sigma
def calculate_sigma(ds):
    ds['sigma'] = gsw.sigma0(ds.vosaline, ds.votemper)
    ds.sigma.attrs['long_name'] = 'Sigma at surface'
    ds.sigma.attrs['units'] = 'kg/m3'
    return ds


# Calculate Practical Salinity, Potential Temperature
def calculate_Sp_theta(ds):
    ds['practical_salinity'] = gsw.conversions.SP_from_SA(ds.vosaline, ds.deptht,
                                                               lon=-126.19, lat=50.48)
    
    ds.practical_salinity.attrs['long_name'] = 'Practical salinity'
    ds.practical_salinity.attrs['units'] = 'None'
    
    ds['potential_temperature'] = gsw.conversions.pt_from_CT(ds.vosaline, ds.votemper)
    
    ds.potential_temperature.attrs['long_name'] = 'Potential temperature'
    ds.potential_temperature.attrs['units'] = 'oC'
    return ds


# Calculate pressure
def __pressure_from_sigma(ds):
    g = 9.81
    pressure = np.zeros_like(ds.deptht.values)
    in_situ_density = ds.sigma[:, :, 0, 0].mean(axis=0) + 1000.
    pressure[0] = (g * in_situ_density[0] * ds.deptht[0]) / 1e4
    for iz in range(1, 40):
        pressure[iz] = pressure[iz-1] + (g * (
            in_situ_density[iz] + in_situ_density[iz-1]) * 0.5 *
            (ds.deptht[iz] - ds.deptht[iz-1])) / 1e4
    return pressure


def  calculate_pressure(ds):
    ds['pressure'] = xr.DataArray(__pressure_from_sigma(ds),
                                  coords={'deptht': ds.deptht},
                                  attrs= {'long_name': 'Pressure',
                                            'units': 'dbar'})
    return ds


def _convert_variable(ds, variable, conversion_uMolar_to_umolkg):
    # Change to umol/kg
    save_attrs = ds[variable].attrs
    ds[variable] = ds[variable] * conversion_uMolar_to_umolkg
    save_attrs['units'] = 'umol/kg'
    ds[variable].attrs = save_attrs
    return ds
    

# Change units from uMol to umol/kg
def change_to_umol_kg(ds):
# Change DIC, TA and OXY to umol/kg from mmol/m3
# any value in mmol/m3 to umol/kg requires
#mmol/m3 * m3/((1000+sigma)kg)*(1000umol/mmol)
    conversion_uMolar_to_umolkg = 1000 / (1000 + ds.sigma)
    if ds.DIC.attrs['units'] != 'umol/kg':
        ds['DIC_uMol'] = ds.DIC
        for variable in ['DIC', 'TA', 'OXY', 'Si', 'NO3']:
            ds =_convert_variable(ds, variable, conversion_uMolar_to_umolkg)
    else:
        print ('Conversion already done, not run again')
    return ds, conversion_uMolar_to_umolkg


# Calculate when isopycnal last at surface
def calculate_year_isopyncal_surface(ds, model_year):
    # parameters for fit
    params0 = 0.1301889490932413
    params1 = 3.8509914822057825
    params2 = 8.301166081413104 + model_year - 2015 

    ds['pycnal_last_at_surface'] = model_year - (params0 * np.exp(-params1 * (25.15 - ds.sigma)) + params2)
    ds.pycnal_last_at_surface.attrs['long_name'] = 'Surface Date'
    ds.pycnal_last_at_surface.attrs['units'] = 'year'
    return ds


# Calculate AOU in Carbon units
def calculate_AOU_in_Carbon(ds):
    # calculate saturated oxygen value, in umol/kg
    ds['osol'] = gsw.O2sol(ds.vosaline, ds.votemper, ds.pressure, lon=-125, lat=50)
    ds.osol.attrs['long_name'] = 'oxygen solubility'
    ds.osol.attrs['units'] = 'umol/kg'
    # calculate AOU in umol/kg
    ds['AOU'] = ds.osol - ds.OXY
    ds.AOU.attrs['long_name'] = 'apparent oxygen utilization'
    ds.AOU.attrs['units'] = 'umol/kg'
    # estimate inorganic carbon  produced by AOU 
    ds['AOU_stoich'] = ds.AOU * (117/170) #Anderson, L. A., & Sarmiento, J. L. (1994)
    ds.AOU_stoich.attrs['long_name'] = 'oxygen solubility in carbon units'
    ds.AOU_stoich.attrs['units'] = 'umol C/kg'
    return ds


# calculate preformed DIC
def calculate_preformed_DIC(ds):
    ds['preformed_DIC'] = ds.DIC - ds.AOU_stoich
    ds.preformed_DIC.attrs['long_name'] = 'preformed dissolved inorganic carbon'
    ds.preformed_DIC.attrs['units'] = 'umol/kg'
    return ds


# calculate preformd pC02
#.   https://pyco2sys.readthedocs.io/en/latest/
#.   and in particular example 2  https://github.com/mvdh7/PyCO2SYS-examples/blob/master/CO2SYSExample2.ipynb
def calculate_preformed_pCO2(ds):
    # Define input and output conditions
    kwargs = dict(
        par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
        par1 = ds.TA.values,  # value of the first parameter
        par2_type = 2,  # The second parameter supplied is of type "2", which means "DIC"
        par2 = ds.preformed_DIC.values,
        salinity = ds.practical_salinity.values,  # Salinity of the sample
        temperature = ds.potential_temperature.values,  # Temperature at input conditions 
        total_silicate = ds.Si.values,  # Concentration of silicate  in the sample (in umol/kg)
        total_phosphate = ds.NO3.values * 1/16.,  # Concentration of phosphate in the sample (in umol/kg)
        opt_k_carbonic = 14,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("14" means Millero 2010")
                              # I checked this in the code
        opt_k_bisulfate = 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
        opt_total_borate = 2,  # Choice of boron:sal ("2" means "LKB10"). # Lee et al what Mucci
        opt_k_fluoride = 1,
        opt_pressured_kCO2 = 1
    )
    # calculate preformed pCO2 - when water was last at surface, ocean pCO2
    results_pCO2 = pyco2.sys(**kwargs)
    ds['preformed_pco2'] = xr.DataArray(results_pCO2['pCO2'],
                                  coords=ds.coords,
                               attrs={'long_name': 'preformed Surface pCO2',
                                        'units': 'uatm'})
    return ds


# calculate disequilibrium pCO2
def calculate_disequ_pCO2(ds):
    ds['diseqPCO2'] = ds.preformed_pco2 - ds.atpco2
    ds.diseqPCO2.attrs['long_name'] = 'ocean - atm, disequilibrium pCO2'
    ds.diseqPCO2.attrs['units'] = 'uatm'
    return ds


# add the PI atmospheric pCO2
def calculate_preformed_PI_pCO2(ds, assumed_pi_pco2):
    ds['PIpref_pco2_inc_diseqpco2'] = ds.diseqPCO2 + assumed_pi_pco2
    ds.PIpref_pco2_inc_diseqpco2.attrs['long_name'] = 'PI preformed pCO2 including disequilbrium'
    ds.PIpref_pco2_inc_diseqpco2.attrs['units'] = 'uatm'
    return ds

   
def calculate_preind_preformed_DIC(ds):
    kwargs = dict(
    par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
    par1 = ds.TA.values,  # value of the first parameter
    par2_type = 4,  # The second parameter supplied is of type "4", which means "pCO2"
    par2 = ds.PIpref_pco2_inc_diseqpco2.values,  # value of the second parameter
    salinity = ds.practical_salinity.values,  # Salinity of the sample
    temperature = ds.potential_temperature.values,  # Temperature at input conditions 
    total_silicate = ds.Si.values,  # Concentration of silicate  in the sample (in umol/kg)
    total_phosphate = ds.NO3.values * 1/16.,  # Concentration of phosphate in the sample (in umol/kg)
    opt_k_carbonic = 14,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("14" means Millero 2010")
    opt_k_bisulfate = 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
    opt_total_borate = 2,  # Choice of boron:sal ("2" means "LKB10")
    opt_k_fluoride = 1,
    opt_pressured_kCO2 = 1
    )
    # Calculate preformed DIC
    dic_results = pyco2.sys(**kwargs)
    ds['preind_pref_dic'] = xr.DataArray(dic_results['dic'],
                                  coords=ds.coords,
                               attrs={'long_name': 'PI preformed DIC',
                                        'units': 'umol/kg'})
    return ds


def calculate_equilibrium_DIC(ds, pCO2, variable, longname, early):
    # only need surface value!
    kwargs = dict(
    par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
    par1 = ds.TA[:, 0].values,  # value of the first parameter
    par2_type = 4,  # The second parameter supplied is of type "4", which means "pCO2"
    par2 = pCO2,  # value of the second parameter
    salinity = ds.practical_salinity[:, 0].values,  # Salinity of the sample
    temperature = ds.potential_temperature[:, 0].values,  # Temperature at input conditions 
    total_silicate = ds.Si[:, 0].values,  # Concentration of silicate  in the sample (in umol/kg)
    total_phosphate = ds.NO3[:, 0].values * 1/16.,  # Concentration of phosphate in the sample (in umol/kg)
    opt_k_carbonic = 14,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("14" means Millero 2010")
    opt_k_bisulfate = 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
    opt_total_borate = 2,  # Choice of boron:sal ("2" means "LKB10")
    opt_k_fluoride = 1,
    opt_pressured_kCO2 = 1
    )
    # Calculate preformed DIC
    dic_results = pyco2.sys(**kwargs)
    if early:
        ds[variable] = xr.DataArray(dic_results['dic'],
                                  coords={'time_counter': ds.time_counter,
                                          'yb': ds.yb,
                                          'xb': ds.xb},
                               attrs={'long_name': longname,
                                        'units': 'umol/kg'})
    else:
        ds[variable] = xr.DataArray(dic_results['dic'],
                                  coords={'time_counter': ds.time_counter,
                                          'yb': ds.yb,
                                          'xbT': ds.xbT},
                               attrs={'long_name': longname,
                                        'units': 'umol/kg'})
    return ds

    
def calculate_delta_DIC(ds):
    ds['deltaDIC'] = ds.preformed_DIC - ds.preind_pref_dic
    ds.deltaDIC.attrs['long_name'] = 'Change in DIC from Preindustrial'
    ds.deltaDIC.attrs['units'] = 'umol/kg'
    return ds


def calculate_final_deep_DIC(ds):
    ds['Final_PI_DIC'] = ds.DIC - ds.deltaDIC
    ds.Final_PI_DIC.attrs['long_name'] = 'Final deep preindustrial DIC'
    ds.Final_PI_DIC.attrs['units'] = 'umol/kg'
    return ds
    

def calculate_dic_diseq(ds):
    # surface only
    ds['dic_diseq'] =  -ds.dic_equ_atmos_now + ds.DIC[:, 0]
    ds.dic_diseq.attrs['long_name'] = 'now DIC difference from equilibrium with now atmos pCO2'
    ds.dic_diseq.attrs['units'] = 'umol/kg'
    return ds


def calculate_surf_DIC_PI(ds):
    ds['surf_DIC_PI'] = ds.dic_equ_atmos_pi + ds.dic_diseq
    ds.surf_DIC_PI.attrs['long_name'] = 'PI Surface DIC'
    ds.surf_DIC_PI.attrs['units'] = 'umol/kg'
    return ds


def correct_near_surface_DIC(ds):
    # Calculate: I don't love this, it has a loop
    surf_intrusion = ds.DIC[0, 0] - ds.surf_DIC_PI  # uM
    sigma_crit = 25 # Value where we switch from surface intrusion to tracked age

    ds['PI_DIC_w_surf_corr'] = xr.DataArray(1.*ds['Final_PI_DIC'].values,
                                  coords=ds.coords,
                               attrs={'long_name': 'PI DIC with Surface Adjustment',
                                        'units': 'umol/kg'})
    depth_crit = np.zeros(ds.sigma.shape[3])
    intrusion_at_crit = np.zeros_like(depth_crit)
    intrusion_values = np.zeros((ds.sigma.shape[1], ds.sigma.shape[3]))

    for ii in range(ds.sigma.shape[3]):
        depth_crit[ii] = np.interp(sigma_crit, ds.sigma[0, :, 0, ii], ds.deptht.values)
        intrusion_at_crit[ii] = np.interp(sigma_crit, ds.sigma[0, :, 0, ii], ds.deltaDIC[0, :, 0, ii])
        intrusion_values[:, ii] = (surf_intrusion[0, ii] * (1 - ds.deptht / depth_crit[ii]) 
                               + intrusion_at_crit[ii] * ds.deptht / depth_crit[ii]);
        kk_crit = np.argmax(ds.deptht.values > depth_crit[ii])
        ds.PI_DIC_w_surf_corr[0, :, 0, ii][ds.deptht < depth_crit[ii]] = (ds.DIC[0, :, 0, ii][ds.deptht < depth_crit[ii]] 
                                                              - intrusion_values[:kk_crit, ii])

    return ds, depth_crit, intrusion_at_crit, intrusion_values


def write_file(ds, date, conversion_uMolar_to_umolkg, early):
    ds['PI_DIC_uMol'] = xr.DataArray(ds.PI_DIC_w_surf_corr.values / conversion_uMolar_to_umolkg,
                                  coords=ds.coords,
                               attrs={'long_name': 'PI DIC with Surface Adjustment in uMol',
                                        'units': 'mmol/m3'})
    new_dataset = xr.DataArray(ds.PI_DIC_uMol, dims=ds.dims,
                                      coords=ds.coords,
                                   attrs=ds.PI_DIC_uMol.attrs).to_dataset()
    ds_attrs = {
        "acknowledgements": "Tereza Jarnikova, Debby Ianson and Susan Allen Algorithms",
        "creator_email": "sallen@eoas.ubc.ca",
        "creator_name": "SalishSeaCast Project Contributors",
        "creator_url": "https://salishsea-meopar-docs.readthedocs.org/",
        "institution": "UBC EOAS",
        "institution_fullname": (
            "Earth, Ocean & Atmospheric Sciences," " University of British Columbia"
        ),
        "summary": (
            "PI DIC: based on age estimate of water, AOU, and near surface delta DIC."
            "NOTE: Date is date of base files, DIC is for circa 1980"
        ),
        "source": (
            "https://github.com/SalishSeaCast/analysis-susan/tree/master/notebooks/Carbon/PI_BOUND_COND.py"
        ),
        "history": (
            "[{}] File creation.".format(datetime.datetime.today().strftime("%Y-%m-%d"))
        ),
    }
    new_dataset.attrs = ds_attrs

    filename = f'./PI_JdF_DIC_v20250708_y{date.year}m{date.month:02d}d{date.day:02d}.nc'
    encoding = {var: {"zlib": True} for var in ['PI_DIC_uMol'] }
    if not early:
        encoding["time_counter"] = {"units": "minutes since 1970-01-01 00:00"}

    new_dataset.to_netcdf(
        path=filename,
        unlimited_dims=("time_counter"),
        encoding=encoding,
    )


