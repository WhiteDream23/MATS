#Import Statements
import numpy as np
import pandas as pd
import os, sys
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")
sns.set_style("ticks")
sns.set_context("poster")
parent_dir = os.path.abspath('../../')
sys.path.append(parent_dir)
import MATS

from MATS.hapi import *
pd.set_option("display.max_rows", 101)


def HITRANlinelist_to_csv(isotopes, minimum_wavenumber, maximum_wavenumber, tablename='tmp', temperature=296,
                          calculate_aw=False):
    #calculate_aw 是否计算理论速度偏移
    """Generates two .csv files generated information available from HTIRAN.  The first line list matches the information available from HITRAN (_HITRAN.csv) and the second supplements the HITRAN information with theoretical values and translates into MATS input format (_initguess.csv)

    Outline

    1. Gets a line list from HITRAN and saves all available parameters to filename_HITRAN.csv
    2. Goes through the data provided from HITRAN and collects the highest order line shape information.
    3.  Where there is missing information for the complete HTP linelist set to 0 or make the following substitutions
        - for missing diluent information fill values with air
        - set missing shift temperature dependences equal to 0 (linear temperature dependence)
        - calculate the SD_gamma based on theory (if calculate aw = True)
        - set the gamma_2 temperature exponent equal to the gamma0 temperature exponent
        - set the delta_2 temperature exponent equal to the delta0 temperature exponent
        - set the dicke narrowing temperature exponent to 1
    4. Save the supplemented and MATS formatted HITRAN information as filename_initguess.csv


    Parameters
    ----------
    isotopes : list
        list of the HITRAN global isotope numbers to include in the HAPI call
    minimum_wavenumber : float
        minimum line center (cm-1) to include in the HAPI call.
    maximum_wavenumber : float
        maximum line center (cm-1) to include in the HAPI call.
    tablename : str, optional
        desired name for table generated from HAPI call. The default is 'tmp'.
    temperature : float, optional
        Nominal temperature of interest.  HITRAN breaks-up the HTP line parameters into temperature regimes.  This allows for selection of the most approriate parameter information. The default is 296.
    calculate_aw : float, optional
        Boolean flag to present option to calculate speed-dependent shift based on theoretical approximation based on temperature exponent, mass of the absorber, and mass of the perturber

    Returns
    -------
    linelist_select : dataframe
        pandas dataframe corresponding to the HITRAN information supplemented by theoretical values/assumptions.
    tablename_HITRAN.csv : .csv file
        file corresponding to available HITRAN information
    tablename_initguess.csv : .csv file
        file corresponding to available HITRAN information supplemented by theory and assumptions in MATS format

    """
    # Possible Parameters, these labels have to match with the HITRANOnline labeling convention. These were taken from HAPI, but this could adapt without updates to HAPI if additional HITRANOnline labels were known.
    HITRAN_parameter_list = ['trans_id', 'molec_id', 'local_iso_id',
                             'nu', 'sw', 'a', 'elower',
                             'gamma_air', 'delta_air', 'n_air', 'deltap_air', 'y_air', 'SD_air', 'Y_SDV_air_296',
                             'beta_g_air',
                             'delta_HT_0_air_296', 'deltap_HT_air_296', 'delta_HT_2_air_296',
                             'nu_HT_air', 'kappa_HT_air', 'eta_HT_air', 'Y_HT_air_296',
                             'gamma_self', 'delta_self', 'n_self', 'deltap_self', 'y_self', 'SD_self', 'Y_SDV_self_296',
                             'beta_g_self',
                             'gamma_HT_0_self_50', 'n_HT_self_50', 'gamma_HT_2_self_50',
                             'delta_HT_0_self_50', 'deltap_HT_self_50', 'delta_HT_2_self_50',
                             'gamma_HT_0_self_150', 'n_HT_self_150', 'gamma_HT_2_self_150',
                             'delta_HT_0_self_150', 'deltap_HT_self_150', 'delta_HT_2_self_150',
                             'gamma_HT_0_self_296', 'n_HT_self_296', 'gamma_HT_2_self_296',
                             'delta_HT_0_self_296', 'deltap_HT_self_296', 'delta_HT_2_self_296',
                             'gamma_HT_0_self_700', 'n_HT_self_700', 'gamma_HT_2_self_700',
                             'delta_HT_0_self_700', 'deltap_HT_self_700', 'delta_HT_2_self_700',
                             'nu_HT_self', 'kappa_HT_self', 'Y_HT_self_296',
                             'gamma_H2', 'delta_H2', 'deltap_H2', 'n_H2',
                             'gamma_CO2', 'delta_CO2', 'n_CO2',
                             'gamma_He', 'delta_He', 'n_He',
                             'gamma_H2O', 'n_H2O',
                             'gp', 'gpp', 'statep', 'statepp', 'global_upper_quanta', 'global_lower_quanta',
                             'local_upper_quanta', 'local_lower_quanta']

    # Retrieves linelist from HITRAN will generate two dataframes, one is the list of everything in HITRAN the other will be an initial guess file for your fitting
    db_begin('data')
    fetch_by_ids(tablename, global_isotopes, minimum_wavenumber, maximum_wavenumber,
                 Parameters=HITRAN_parameter_list)  # pulls down all HITRAN data for the data #, ParameterGroups=PARLIST_ALL, , ParameterGroups = [PARAMETER_GROUPS['standard']]
    cond = ('AND', ('between', 'nu', minimum_wavenumber, maximum_wavenumber), ('>=', 'sw', intensity_cutoff))
    select(tablename, Conditions=cond, DestinationTableName='tmp')

    # Generates Pandas dataframe using linelist
    ## First Generates a table that contains everything that is in HITRAN for that molecule in that waverange
    linelist = pd.DataFrame()
    standard = ['trans_id', 'molec_id', 'local_iso_id', 'nu', 'sw', 'a', 'elower', 'gp', 'gpp', 'elower',
                'global_upper_quanta', 'global_lower_quanta', 'local_upper_quanta', 'local_lower_quanta']
    for par_name in standard:  # adds standard items to the linelist
        linelist[par_name] = LOCAL_TABLE_CACHE['tmp']['data'][par_name.lower()]
    for par_name in PARLIST_HT_AIR:
        try:
            linelist[par_name] = LOCAL_TABLE_CACHE['tmp']['data'][par_name.lower()]
        except:
            pass
    for par_name in PARLIST_HT_SELF:
        try:
            linelist[par_name] = LOCAL_TABLE_CACHE['tmp']['data'][par_name.lower()]
        except:
            pass

    for par_name in PARLIST_VOIGT_AIR:  # adds VP air parameters
        try:
            linelist[par_name] = LOCAL_TABLE_CACHE['tmp']['data'][par_name.lower()]
        except:
            pass

    for par_name in PARLIST_VOIGT_SELF:  # adds VP self parameters
        try:
            linelist[par_name] = LOCAL_TABLE_CACHE['tmp']['data'][par_name.lower()]
        except:
            pass

    for par_name in PARLIST_VOIGT_H2:
        try:
            linelist[par_name] = LOCAL_TABLE_CACHE['tmp']['data'][par_name.lower()]
        except:
            pass
    for par_name in PARLIST_VOIGT_CO2:
        try:
            linelist[par_name] = LOCAL_TABLE_CACHE['tmp']['data'][par_name.lower()]
        except:
            pass
    for par_name in PARLIST_VOIGT_HE:
        try:
            linelist[par_name] = LOCAL_TABLE_CACHE['tmp']['data'][par_name.lower()]
        except:
            pass
    for par_name in PARLIST_VOIGT_H2O:
        try:
            linelist[par_name] = LOCAL_TABLE_CACHE['tmp']['data'][par_name.lower()]
        except:
            pass
    for par_name in ['y_air', 'y_self']:
        try:
            linelist[par_name] = LOCAL_TABLE_CACHE['tmp']['data'][par_name.lower()]
        except:
            pass
    for par_name in ['SD_air', 'SD_self', 'SD_H2', 'SD_CO2', 'SD_He',
                     'SD_H2O']:  # Addition of Speed Dependent Parameters
        try:
            linelist[par_name] = LOCAL_TABLE_CACHE['tmp']['data'][par_name.lower()]
        except:
            pass
    for par_name in ['beta_g_air', 'beta_g_self', 'beta_g_H2', 'beta_g_CO2', 'beta_g_He', 'beta_g_H2O']:
        try:
            linelist[par_name] = LOCAL_TABLE_CACHE['tmp']['data'][par_name.lower()]
        except:
            pass

    for par_name in list(linelist):
        if all(isinstance(val, np.ma.core.MaskedConstant) for val in linelist[par_name]):
            linelist.drop(columns=par_name, inplace=True)

    linelist.to_csv(tablename + '_HITRAN.csv', index=False)
    # Next segment looks at at all possible data and makes an initial guess fill prioritizing HTP over non HTP parameters, but will mix them.
    # It will also calculate aw based on theory and set values where there is an air value, but no self value equal to the air value.  COmment line printed at end will detail these
    #  Will also use GP values if available (and NGP from HTP are not) will make a note in comment if this happens
    avail_species = []
    linelist_select = linelist[
        ['trans_id', 'molec_id', 'local_iso_id', 'nu', 'sw', 'a', 'elower', 'gp', 'gpp', 'elower',
         'global_upper_quanta', 'global_lower_quanta', 'local_upper_quanta', 'local_lower_quanta']].copy()
    for param in list(linelist): # 1105 air\self
        if ('gamma_' in param) and ('HT' not in param):
            avail_species.append(param[6:])

    # define reference temperature and pressure
    Tref = 296.  # K

    # define actual temperature and pressure
    T = temperature  # K
    TRanges = [(0, 100), (100, 200), (200, 400), (400, float('inf'))]
    Trefs = [50., 150., 296., 700.]
    for TRange, TrefHT in zip(TRanges, Trefs):
        if T >= TRange[0] and T < TRange[1]:
            break

    for molecule in linelist_select['molec_id'].unique():
        for isotope in linelist_select[linelist_select['molec_id'] == molecule]['local_iso_id'].unique():
            for species in avail_species:
                comment = ''
                # Gamma0
                #使用HITRANOnline中的gamma_air和gamma_self
                try:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'gamma0_%s' % (species)] = \
                    linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                        'gamma_HT_0_%s_%d' % (species, TrefHT)].values
                except:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'gamma0_%s' % (species)] = \
                    linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                        'gamma_%s' % species].values
                # Temperature Dependence Gamma0
                #使用HITRANOnline中的n_air和n_self
                try:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'n_gamma0_%s' % (species)] = \
                    linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                        'n_HT_%s_%d' % (species, TrefHT)].values
                except:
                    try:
                        linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                    linelist_select['local_iso_id'] == isotope), 'n_gamma0_%s' % (species)] = \
                        linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                            'n_%s' % species].values
                    except:
                        linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                    linelist_select['local_iso_id'] == isotope), 'n_gamma0_%s' % (species)] = \
                        linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                            'n_air'].values
                if (linelist_select[
                        (linelist_select['molec_id'] == molecule) & (linelist_select['local_iso_id'] == isotope)][
                        'n_gamma0_%s' % (species)] == 0).all():
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'n_gamma0_%s' % (species)] = \
                    linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)]['n_air'].values
                    comment += 'set n_gamma0_%s' % (species) + ' to n_gamma0_air'

                # Delta0
                try:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'delta0_%s' % (species)] = \
                    linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                        'delta_HT_0_%s_%d' % (species, TrefHT)].values
                except:
                    try:
                        linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                    linelist_select['local_iso_id'] == isotope), 'delta0_%s' % (species)] = \
                        linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                            'delta_%s' % species].values
                    except:
                        linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                    linelist_select['local_iso_id'] == isotope), 'delta0_%s' % (species)] = \
                        linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                            'delta_air']
                        comment += 'set delta0_%s' % (species) + ' to delta0_air'
                # Temperature Dependence of Delta0
                try:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'n_delta0_%s' % (species)] = \
                    linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                        'deltap_HT_%s_%d' % (species, TrefHT)].values
                except:
                    try:
                        linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                    linelist_select['local_iso_id'] == isotope), 'n_delta0_%s' % (species)] = \
                        linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                            'deltap_%s' % species].values
                    except:
                        try:
                            linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                        linelist_select['local_iso_id'] == isotope), 'n_delta0_%s' % (species)] = \
                            linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                                'deltap_air'].values
                            comment += 'set n_delta0_%s' % (species) + ' to n_delta0_air'
                        except:
                            linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                        linelist_select['local_iso_id'] == isotope), 'n_delta0_%s' % (species)] = 0

                # Speed Dependent Broadening
                try:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'SD_gamma_%s' % (species)] = \
                    linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                        'gamma_HT_2_%s_%d' % (species, TrefHT)].values
                except:
                    try:
                        linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                    linelist_select['local_iso_id'] == isotope), 'SD_gamma_%s' % (species)] = \
                        linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                            'SD_%s' % species].values
                    except:
                        try:
                            linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                        linelist_select['local_iso_id'] == isotope), 'SD_gamma_%s' % (species)] = \
                            linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                                'SD_air'].values
                        except:

                            if calculate_aw:
                                for i in linelist_select[(linelist_select['molec_id'] == molecule) & (
                                        linelist_select['local_iso_id'] == isotope)].index:
                                    if species == 'air':
                                        m_p = 28.97
                                    elif species == 'H2':
                                        m_p = 2.01588
                                    elif species == 'CO2':
                                        m_p = 43.98983
                                    elif species == 'HE':
                                        m_p = 4.002602
                                    elif species == 'H2O':
                                        m_p = 18.010565
                                    else:
                                        m_p = float(
                                            molecularMass(linelist.loc[i]['molec_id'], linelist.loc[i]['local_iso_id']))
                                    m_a = float(
                                        molecularMass(linelist.loc[i]['molec_id'], linelist.loc[i]['local_iso_id']))
                                    aw = (1 - linelist_select.loc[i, 'n_gamma0_%s' % (species)]) * (2 / 3) * (
                                                (m_p / m_a) / (1 + (m_p / m_a)))
                                    linelist_select.loc[i, 'SD_gamma_%s' % (species)] = aw
                            else:
                                linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                            linelist_select['local_iso_id'] == isotope), 'SD_gamma_%s' % (species)] = 0

                # Temperature Dependence of SD Broadening
                linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                            linelist_select['local_iso_id'] == isotope), 'n_gamma2_%s' % (species)] = linelist_select[
                    (linelist_select['molec_id'] == molecule) & (linelist_select['local_iso_id'] == isotope)][
                    'n_gamma0_%s' % (species)].values

                # Speed Dependent Shift
                try:
                    linelist_select.loc[linelist_select.index(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'SD_delta_%s' % (species)] = \
                    linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                        'delta_HT_2_%s_%d' % (species, TrefHT)].values
                except:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'SD_delta_%s' % (species)] = 0
                # Temperature Dependence of SD Shift
                linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                            linelist_select['local_iso_id'] == isotope), 'n_delta2_%s' % (species)] = \
                linelist_select[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                    'n_delta0_%s' % (species)].values

                # nuVC
                try:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'nuVC_%s' % (species)] = \
                    linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                        'nu_HT_%s' % species].values
                except:
                    try:
                        linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                    linelist_select['local_iso_id'] == isotope), 'nuVC_%s' % (species)] = \
                        linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                            'beta_g_%s' % species].values
                        comment += ' nuVC_' + species + ' = beta_g_' + species
                    except:
                        linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                    linelist_select['local_iso_id'] == isotope), 'nuVC_%s' % (species)] = 0
                # nuVC temperature dependence
                try:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'n_nuVC_%s' % (species)] = \
                    linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                        'kappa_HT_%s' % species].values
                except:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'n_nuVC_%s' % (species)] = 1

                # eta
                try:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'eta_%s' % (species)] = \
                    linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                        'eta_HT_%s' % species].values
                except:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'eta_%s' % (species)] = 0

                    # Linemixing
                try:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'y_%s' % (species) + '_' + str(
                        int(temperature))] = \
                    linelist[(linelist['molec_id'] == molecule) & (linelist['local_iso_id'] == isotope)][
                        'y_%s' % species].values
                except:
                    linelist_select.loc[(linelist_select['molec_id'] == molecule) & (
                                linelist_select['local_iso_id'] == isotope), 'y_%s' % (species) + '_' + str(
                        int(temperature))] = 0

                print(molecule, isotope, species, comment)
    linelist_select.to_csv(tablename + '_initguess.csv', index=False)

    return linelist_select

tablename = 'CO2'
global_isotopes = [7,8,9,10]
wave_min = 9505
wave_max = 9507
intensity_cutoff = 1e-30

linelist_select = (HITRANlinelist_to_csv(global_isotopes, wave_min, wave_max, tablename = tablename, calculate_aw = True))