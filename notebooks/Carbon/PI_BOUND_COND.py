import gsw
import numpy as np
import pandas as pd
import PyCO2SYS as pyco2


# two functions to help calculate tco2
def _co2_from_year_init():
    co2_rec = pd.read_csv('../../../analysis_tereza/notebooks/PI_CARBON_PAPER/PI_BOUND_COND/CLEAN/lawdome_maunaloa.csv')
    sindex = co2_rec[co2_rec['YEAR'] == 1950].index.values[0]
    co2_rec = co2_rec.set_index('YEAR')
    c = np.polynomial.Polynomial.fit(co2_rec.index[sindex:], co2_rec.PPMCO2[sindex:], deg=2)
    return c
    

def _co2_from_year(ds, c):
    if (np.min(ds.pycnal_last_at_surface.values) < 1950) | (np.max(ds.pycnal_last_at_surface.values) > 2025):
        print ('Years asked for are out of range of fit!', np.min(ds.pycnal_last_at_surface.values),
              np.max(ds.pycnal_last_at_surface.values))
    ds = ds.assign(atpco2=ds.pycnal_last_at_surface)
    ds.atpco2[:] = c(ds.pycnal_last_at_surface.values)
    return ds
    

# Calculate Sigma
def calculate_sigma(ds):
    ds['sigma'] = gsw.sigma0(ds.vosaline, ds.votemper)
    ds.sigma.attrs['long_name'] = 'Sigma at surface'
    ds.sigma.attrs['units'] = 'kg/m3'
    return ds


def _convert_variable(ds, variable):
    # Change to umol/kg
    conversion_uMolar_to_umolkg = 1000 / (1000 + ds.sigma)
    save_attrs = ds.[variable].attrs
    ds[variable] = ds[variable] * conversion_uMolar_to_umolkg
    save_attrs['units'] = 'umol/kg'
    ds[variable].attrs = save_attrs
    return ds
    

# Change units from uMol to umol/kg
def change_to_umol_kg(ds):
# Change DIC, TA and OXY to umol/kg from mmol/m3
# any value in mmol/m3 to umol/kg requires
#mmol/m3 * m3/((1000+sigma)kg)*(1000umol/mmol)
    if ds.DIC.attrs['units'] != 'umol/kg':
        ds['DIC_uMol'] = ds.DIC
        for variable in ['DIC', 'TA', 'OXY', 'Si', 'NO3']:
            ds=_convert(ds, variable,conversion_uMolar_to_umolkg )
    else:
        print ('Conversion already done, not run again')
    return ds


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


# # Calculate the Atmospheric CO2 when watermass was at surface: atpco2
def calculate_atpco2_for_year(ds):
    c = _co2_from_year_init()
    ds = _co2_from_year(ds, c)  # uses the Mauna Loa look up
    ds.atpco2.attrs['long_name'] = 'Atmospheric pCO2'
    ds.atpco2.attrs['units'] = 'ppm'
    return ds


# Calculate AOU in Carbon units
def calculate_AOU_in_Carbon(ds):
    # calculate saturated oxygen value, pressure is ds.deptht
    ds['osol'] = gsw.O2sol(ds.vosaline, ds.votemper, ds.deptht, lon=-125, lat=50)
    ds.osol.attrs['long_name'] = 'oxygen solubility'
    ds.osol.attrs['units'] = 'umol/kg'
    # calculate AOU in umol/kg
    ds['AOU'] = ds.osol - ds. OXY
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
def calculate_preformed_pCO2(ds, conversion_uMolar_to_umolkg):
    # Translate Salinity and Conservative Temperature to psu and potential temperature
    ds['psu'] = gsw.SP_from_SA(ds.vosaline, p=0*ds.vosaline, lat=49, lon=-126)
    ds['pottemp'] = gsw.pt_from_CT(ds.vosaline, ds.votemper)
    # Define input and output conditions
    kwargs = dict(
        par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
        par1 = ds.TA.values,  # value of the first parameter
        par2_type = 2,  # The second parameter supplied is of type "2", which means "DIC"
        par2 = ds.preformed_DIC.values,
        salinity = ds.psu.values,  # Salinity of the sample
        temperature = ds.pottemp.values,  # Temperature at input conditions 
        total_silicate = 0*ds.Si.values * conversion_uMolar_to_umolkg,  # Concentration of silicate  in the sample (in umol/kg)
        total_phosphate = 0*ds.NO3.values * conversion_uMolar_to_umolkg * 1/16.,  # Concentration of phosphate in the sample (in umol/kg)
        opt_k_carbonic = 14,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("14" means Millero 2010")
                              # I checked this in the code
        opt_k_bisulfate = 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
        opt_total_borate = 2,  # Choice of boron:sal ("2" means "LKB10"). # Lee et al what Mucci
        opt_k_fluoride = 1,
        opt_pressured_kCO2 = 1
    )
    # calculate preformed pCO2 - when water was last at surface, ocean pCO2
    results_pCO2 = pyco2.sys(**kwargs)
    ds = ds.assign(preformed_pco2= 0.*ds['atpco2'])
    ds.preformed_pco2[:] = results_pCO2['pCO2']
    ds.preformed_pco2.attrs['long_name'] = 'preformed Surface pCO2'
    ds.preformed_pco2.attrs['units'] = 'uatm'    
    return ds


# calculate disequilibrium pCO2
def calculate_disequ_pCO2(ds):
    ds['diseqPCO2'] = ds.preformed_pco2 - ds.atpco2
    ds.diseqPCO2.attrs['long_name'] = 'ocean - atm, disequilibrium pCO2'
    ds.diseqPCO2.attrs['units'] = 'uatm'
    return ds


# add the PI atmospheric pCO2
def calculate_preformed_PI_pCO2(ds):
    ds['PIpref_pco2_inc_diseqpco2'] = ds.diseqPCO2 + 284
    ds.PIpref_pco2_inc_diseqpco2.attrs['long_name'] = 'PI preformed pCO2 including disequilbrium'
    ds.PIpref_pco2_inc_diseqpco2.attrs['units'] = 'uatm'
    return ds


    
def calculate_preformed_What(ds, conversion_uMolar_to_umolkg):
    kwargs = dict(
    par1_type = 1,  # The first parameter supplied is of type "1", which means "alkalinity"
    par1 = ds.TA.values,  # value of the first parameter
    par2_type = 4,  # The second parameter supplied is of type "4", which means "pCO2"
    par2 = ds.PIpref_pco2_inc_diseqpco2.values,  # value of the second parameter
    salinity = ds.psu.values,  # Salinity of the sample
    temperature = ds.pottemp.values,  # Temperature at input conditions 
    total_silicate = 0*ds.Si.values * conversion_uMolar_to_umolkg,  # Concentration of silicate  in the sample (in umol/kg)
    total_phosphate = 0*ds.NO3.values * conversion_uMolar_to_umolkg * 1/16.,  # Concentration of phosphate in the sample (in umol/kg)
    opt_k_carbonic = 14,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("14" means Millero 2010")
    opt_k_bisulfate = 1,  # Choice of HSO4- dissociation constant KSO4 ("1" means "Dickson")
    opt_total_borate = 2,  # Choice of boron:sal ("2" means "LKB10")
    opt_k_fluoride = 1,
    opt_pressured_kCO2 = 1
    )
    # Calculate preformed DIC
    dic_results = pyco2.sys(**kwargs)
    
    ds = ds.assign(preind_pref_dic= 0.*ds['atpco2'])
    ds.preind_pref_dic[:] = dic_results['dic']
    ds.preind_pref_dic.attrs['long_name'] = 'PI preformed DIC'
    ds.preind_pref_dic.attrs['units'] = 'umol/kg'



