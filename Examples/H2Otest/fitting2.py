import numpy as np
import pandas as pd
import os, sys
import matplotlib.pyplot as plt
from matplotlib import gridspec

from MATS import linelistdata

#from MATS.linelistdata import linelistdata
parent_dir = os.path.abspath('../../')
sys.path.append(parent_dir)
import MATS

import seaborn as sns
sns.set_style("whitegrid")
sns.set_style("ticks")
sns.set_context("poster")


#Generic Fit Parameters
wave_range = 1.5 #range outside of experimental x-range to simulate
IntensityThreshold = 1e-30 #intensities must be above this value to be simulated
Fit_Intensity = 1e-29 #intensities must be above this value for the line to be fit
order_baseline_fit = 1
tau_column = 'aaaa' # Mean tau/ us
freq_column = 'Total Frequency /MHz' # Total Frequency /MHz
pressure_column = 'Cavity Pressure /Torr'
temperature_column = 'Cavity Temperature Side 2 /C'


#Define all Spectra individually
spec_1 = MATS.Spectrum('H2O-line1',
                        molefraction = { 1 :1}, natural_abundance = True, diluent = 'self',
                        etalons = {}, baseline_order = 1,
                        input_freq = False, frequency_column = freq_column,
                        input_tau = False, tau_column = tau_column,
                        pressure_column = pressure_column, temperature_column = temperature_column,
                        nominal_temperature = 296, x_shift = 0.00)
spec_2 = MATS.Spectrum('H2O-line2',
                        molefraction = { 1 :1}, natural_abundance = True, diluent = 'self',
                        etalons = {}, baseline_order = 1,
                        input_freq = False, frequency_column = freq_column,
                        input_tau = False, tau_column = tau_column,
                        pressure_column = pressure_column, temperature_column = temperature_column,
                        nominal_temperature = 296, x_shift = 0.00)
spec_3 = MATS.Spectrum('H2O-line3',
                        molefraction = { 1 :1}, natural_abundance = True, diluent = 'self',
                        etalons = {}, baseline_order = 1,
                        input_freq = False, frequency_column = freq_column,
                        input_tau = False, tau_column = tau_column,
                        pressure_column = pressure_column, temperature_column = temperature_column,
                        nominal_temperature = 296, x_shift = 0.00)
spec_4 = MATS.Spectrum('H2O-line4',
                        molefraction = { 1 :1}, natural_abundance = True, diluent = 'self',
                        etalons = {}, baseline_order = 1,
                        input_freq = False, frequency_column = freq_column,
                        input_tau = False, tau_column = tau_column,
                        pressure_column = pressure_column, temperature_column = temperature_column,
                        nominal_temperature = 296, x_shift = 0.00)
#spec_1.plot_wave_alpha()
#Read in linelists
PARAM_LINELIST = linelistdata['H2O_initguess']
#Add all spectrum to a Dataset object
#做了格式统一的处理
SPECTRA = MATS.Dataset([spec_1,spec_2,spec_3,spec_4], 'Line Intensity',PARAM_LINELIST)

#Generate Baseline Parameter list based on number of etalons in spectra definitions and baseline order
BASE_LINELIST = SPECTRA.generate_baseline_paramlist()

FITPARAMS = MATS.Generate_FitParam_File(SPECTRA, PARAM_LINELIST, BASE_LINELIST, lineprofile = 'VP', linemixing = True,
                                  fit_intensity = Fit_Intensity, threshold_intensity = IntensityThreshold, sim_window = wave_range,
                                  nu_constrain = True, sw_constrain = True, gamma0_constrain = True, delta0_constrain = True,
                                   aw_constrain = True, as_constrain = True,
                                   nuVC_constrain = True, eta_constrain =True, linemixing_constrain = True,
                                        )
                                     #additional_columns = ['trans_id', 'local_lower_quanta'])

FITPARAMS.generate_fit_param_linelist_from_linelist(vary_nu = {1:{1:True}}, vary_sw = {1:{1:True}},
                                                    vary_gamma0 = {1:{1:True}}, vary_n_gamma0 = {},
                                                    vary_delta0 = {1:{1:True}}, vary_n_delta0 = {},
                                                    vary_aw = {}, vary_n_gamma2 = {},
                                                    vary_as = {}, vary_n_delta2 = {},
                                                    vary_nuVC = {}, vary_n_nuVC = {},
                                                    vary_eta = {}, vary_linemixing = {1:{1:True}})

FITPARAMS.generate_fit_baseline_linelist(vary_baseline = True, vary_molefraction = {1:False}, vary_xshift = False,
                                      vary_etalon_amp= False, vary_etalon_period= False, vary_etalon_phase= False)
fit_data = MATS.Fit_DataSet(SPECTRA, 'Baseline_LineList', 'Parameter_LineList',
                            minimum_parameter_fit_intensity=Fit_Intensity, weight_spectra=False,
                            baseline_limit=False, baseline_limit_factor=10,
                            molefraction_limit=False, molefraction_limit_factor=1.1,
                            etalon_limit=False, etalon_limit_factor=2,  # phase is constrained to +/- 2pi,
                            x_shift_limit=False, x_shift_limit_magnitude=0.5,
                            nu_limit=True, nu_limit_magnitude=0.1,
                            sw_limit=True, sw_limit_factor=2,
                            gamma0_limit=True, gamma0_limit_factor=3, n_gamma0_limit=False, n_gamma0_limit_factor=50,
                            delta0_limit=True, delta0_limit_factor=2, n_delta0_limit=False, n_delta0_limit_factor=50,
                            SD_gamma_limit=False, SD_gamma_limit_factor=2, n_gamma2_limit=False,
                            n_gamma2_limit_factor=50,
                            SD_delta_limit=False, SD_delta_limit_factor=50, n_delta2_limit=False,
                            n_delta2_limit_factor=50,
                            nuVC_limit=False, nuVC_limit_factor=2, n_nuVC_limit=False, n_nuVC_limit_factor=50,
                            eta_limit=False, eta_limit_factor=50, linemixing_limit=True, linemixing_limit_factor=50)
params = fit_data.generate_params()

#两次打印生成的出来的参数格式是一样的，部分参数的值不一样，fit_data做了更新拟合参数的事
#residual_analysis中绘制了残差图，indv_resid_plot=True表示绘制每个光谱的残差图
#update_params更新了拟合参数，更新了两文件
print(params.pretty_print())
result = fit_data.fit_data(params)
print(result.params.pretty_print())
print(result.success)
print(result.message)
print(result.covar)
fit_data.residual_analysis(result, indv_resid_plot= False)
fit_data.update_params(result,'base','para')
SPECTRA.generate_summary_file(save_file=True)
#SPECTRA.plot_model_residuals()

