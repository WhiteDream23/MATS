#Import Packages
import numpy as np
import pandas as pd
import scipy.interpolate as interpolate
from pathlib import Path

def o2_cia_karman_model(wavenumbers, T, P,Diluent,
                        SO_O2, SO_N2, EXCH_O2, 
                        EXCH_b, EXCH_c, 
                        SO_b_O2_O2, SO_c_O2_O2,
                        SO_b_O2_N2, SO_c_O2_N2,
                        SO_shift_O2_O2 = 0, SO_shift_O2_N2 = 0,
                        EXCH_shift = 0,
                        band = 'singlet_delta'):
    '''
    #Default values based on Karman, T. et al., Icarus 2019, 328 , 160 175.
    
    #Singlet Delta
    #EXCH_c, EXCH_b, EXCH_a = [3.6307626466573398e-06, 0.0028385240774561797, 1]
    #SO_c, SO_b, SO_a =[1.4670403122287775e-06, 0.00014594154382655564, 1]
    #SO_O2, SO_N2, EXCH_O2 = [39.13, 70.74 ,304.7448171031378] #Initial guess derived based HITRAN 2020 reported theoretical CIA
    
    #ABand
    #EXCH_c, EXCH_b, EXCH_a =[6.559060698261758e-05, 0.011869752199984616, 1]
    #SO_c, SO_b, SO_a =[1.5906417750834962e-06, 0.00011263534228667677, 1]
    #SO_O2, SO_N2, EXCH_O2 = [6.20731994222978,7.961801018674746, 39.42079598436756]
    '''
    
    
    prefix_local = Path(__file__).parent / "CIA_Data"
    paths_default = list(prefix_local.glob('*.csv'))
    paths_dict = {p.with_suffix('').name : p for p in paths_default}
    paths_dict['singlet_delta'] = paths_dict.pop('Karman_SingletDelta_Mech_Temp_Dep')
    paths_dict['a_band'] = paths_dict.pop('Karman_ABand_Mech_Temp_Dep')
    
    
    if band == 'singlet_delta':
        df_Karman = pd.read_csv(paths_dict['singlet_delta'])
    elif band == 'a_band':
        df_Karman = pd.read_csv(paths_dict['a_band'])
    #print (len(df_Karman['Wavenumber (cm-1)'].values), len(df_Karman['SO Temp Dep - d'].values))
    
    model_wavenumber = df_Karman['Wavenumber (cm-1)'].values
    SO = (df_Karman['SO Temp Dep - d'].values*(T-296)**3 +df_Karman['SO Temp Dep - c'].values*(T-296)**2 + df_Karman['SO Temp Dep - b'].values*(T-296) + df_Karman['SO Temp Dep - a'].values)*df_Karman['Spin Orbit (cm-1 amagat-2) Normalized'].values    
    EXCH = (df_Karman['EXCH Temp Dep - d'].values*(T-296)**3 +df_Karman['EXCH Temp Dep - c'].values*(T-296)**2 + df_Karman['EXCH Temp Dep - b'].values*(T-296) + df_Karman['EXCH Temp Dep - a'].values)*df_Karman['Exchange (cm-1 amagat-2) Normalized'].values

    
    f = interpolate.interp1d(model_wavenumber, SO, fill_value = 'extrapolate', bounds_error = False, kind = 'slinear')
    SO_O2_O2 = f(wavenumbers- SO_shift_O2_O2)
    g = interpolate.interp1d(model_wavenumber, SO, fill_value = 'extrapolate', bounds_error = False, kind = 'slinear')
    SO_O2_N2 = g(wavenumbers- SO_shift_O2_N2)
    h = interpolate.interp1d(model_wavenumber, EXCH, fill_value = 'extrapolate', bounds_error = False, kind = 'slinear')
    EXCH_ = h(wavenumbers - EXCH_shift)
    
    model_O2_O2 = EXCH_O2*(1 + EXCH_b*(T-296) + EXCH_c*(T-296)**2)*EXCH_ + SO_O2*(1 + SO_b_O2_O2*(T-296) + SO_c_O2_O2*(T-296)**2)*SO_O2_O2
    model_O2_N2 = SO_N2*(1 + SO_b_O2_N2*(T-296) + SO_c_O2_N2*(T-296)**2)*SO_O2_N2
    
    amagats_O2 = (P/1)*(273.15/(T))*Diluent['O2']['composition']
    amagats_N2 = (P/1)*(273.15/(T))*Diluent['N2']['composition']
    
    CIA_model = (model_O2_O2*amagats_O2*amagats_O2) + (model_O2_N2*amagats_N2*amagats_O2) 
    #CIA_model = (Diluent['O2']['composition']*model_O2_O2) + (Diluent['N2']['composition']*model_O2_N2) # in amagats
    
    return CIA_model

