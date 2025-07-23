import arrow
import datetime as dt
import pandas as pd
from pathlib import Path
import xarray as xr

working_dir = Path('/home/sallen/MEOPAR/Karyn_Reshapr/')

CSis, CSie, CSjs, CSje = 450, 500, 250, 300
JFis, JFie, JFjs, JFje = 300, 365, 50, 100
Nks, Nke = 0, 10
Pks, Pke = 0, 24
Zks, Zke = 0, 40


def set_geometry():
    mesh = xr.open_dataset('/home/sallen/MEOPAR/grid/mesh_mask202108.nc')
    mesh = mesh.rename_dims({'z': 'depth', 'x': 'gridX', 'y': 'gridY'})
    # Nitrate for Central SoG
    e3t_CSOG = mesh.e3t_0[0, Nks:Nke, CSis:CSie, CSjs:CSje]
    tmask = mesh.tmask[0, Nks:Nke, CSis:CSie, CSjs:CSje]
    # Nitrate for JdF
    e3t_JF = mesh.e3t_0[0, Nks:Nke, JFis:JFie, JFjs:JFje]
    tmask_JF = mesh.tmask[0, Nks:Nke, JFis:JFie, JFjs:JFje]
    # Phytoplankton for CSOG
    e3t_CSOG_50 = mesh.e3t_0[0, Pks:Pke, CSis:CSie, CSjs:CSje]
    tmask_CSOG_50 = mesh.tmask[0, Pks:Pke, CSis:CSie, CSjs:CSje]
    # Phytoplankton for JdF
    e3t_JF_50 = mesh.e3t_0[0, Pks:Pke, JFis:JFie, JFjs:JFje]
    tmask_JF_50 = mesh.tmask[0, Pks:Pke, JFis:JFie, JFjs:JFje]
    # Zooplankton for CSOG
    e3t_CSOG_AD = mesh.e3t_0[0, :, CSis:CSie, CSjs:CSje]
    tmask_CSOG_AD = mesh.tmask[0, :, CSis:CSie, CSjs:CSje]
    # Zooplankton for JdF
    e3t_JF_AD = mesh.e3t_0[0, :, JFis:JFie, JFjs:JFje]
    tmask_JF_AD = mesh.tmask[0, :, JFis:JFie, JFjs:JFje]

    return (e3t_CSOG, tmask, e3t_JF, tmask_JF, e3t_CSOG_50, tmask_CSOG_50, e3t_JF_50, tmask_JF_50, 
            e3t_CSOG_AD, tmask_CSOG_AD, e3t_JF_AD, tmask_JF_AD)


def find_variances (stem, variables, variable_names, e3t, tmask):
    start = dt.datetime(2007, 1, 1)
    end = dt.datetime(2022, 12, 31)

    columns = ['Year', 'Month'] + variable_names 
    for variable_name in variable_names:
        columns = columns + [f'{variable_name} time', f'{variable_name} space']
    print (columns)
    variance = pd.DataFrame(None, columns=columns)

    for imonth, run_s in enumerate(arrow.Arrow.range('month', start, end)):
        file_path = working_dir.joinpath(
            f'{stem}_{run_s.format("YYYYMMDD")}_{run_s.shift(months=+1, days=-1).format("YYYYMMDD")}.nc')
        print (file_path)
        data = xr.open_dataset(file_path)
        new_row_dict = {'Year': run_s.year, 'Month': run_s.month}
        for variable, variable_name in zip(variables, variable_names):
            new_row_dict[variable_name] = ((((data * e3t).where(tmask ==1).sum(dim='depth'))
                                           /e3t.where(tmask==1).sum(dim='depth')
                                            ).var()[variable].item())     
            new_row_dict[f'{variable_name} time'] = ((((data * e3t).where(tmask ==1).sum(dim='depth'))
                                           /e3t.where(tmask==1).sum(dim='depth')
                                            ).mean(dim='gridX').mean(dim='gridY').var(dim='time')[variable].item())              
            new_row_dict[f'{variable_name} space'] = ((((data * e3t).where(tmask ==1).sum(dim='depth'))
                                           /e3t.where(tmask==1).sum(dim='depth')
                                            ).mean(dim='time').var()[variable].item())              
        variance.loc[imonth] = new_row_dict
        data.close()
    
    return variance


def run_a_case():
    (e3t_CSOG, tmask, e3t_JF, tmask_JF, e3t_CSOG_50, tmask_CSOG_50, e3t_JF_50, tmask_JF_50, 
     e3t_CSOG_AD, tmask_CSOG_AD, e3t_JF_AD, tmask_JF_AD) = set_geometry()

    stem = 'JdF_zoo'
    variable = ['microzooplankton', 'mesozooplankton']
 
    variance = find_variances(stem, variable, ['microzooplankton variance', 'mesozooplankton variance'], e3t_JF_AD, tmask_JF_AD)

    variance.to_csv('JdF_zoo_2007_2022.csv')


if __name__ == "__main__":
    print("Let's Go")
    run_a_case()

    

    