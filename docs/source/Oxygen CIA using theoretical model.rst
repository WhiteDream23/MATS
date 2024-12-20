Oxygen CIA using theoretical model
==================================

Provided in MATS are several examples highlighting MATS capabilities, which can be found in the MATS `examples folder <https://github.com/usnistgov/MATS/tree/master/Examples>`_. 

Karman et al published a theoretical O_2 - O_2 and O_2 - N_2 collision-induced absorption model based on quantum calculations https://doi.org/10.1038/s41557-018-0015-x. This model is composed of spin-orbit and exchange mechanisms. In the manuscript Parameterized Model to Approximate Theoretical Collision-Induced Absorption Band Shapes for O2-O2 and O2-N2, Adkins et al present a parameterized representation of the theoretical model reported by Karman et al.

The details of the model can be found in https://doi.org/10.1016/j.jqsrt.2023.108732.

We have included this parameterized model in MATS to allow for use in multi-spectrum fitting analyses in order to fit monomer and CIA components to the absorption.  `The example below shows a use example for this model. <https://github.com/usnistgov/MATS/tree/master/Examples/O2CIA_Karman_etal_model>`_ 

 Simulating Spectrum
 +++++++++++++++++++
 
 Most of this example follows from the theoretical spectra generation and fit example. However, this example makes use of the new parameterized CIA model based on the theory published by Karman et al https://doi.org/10.1016/j.jqsrt.2023.108732.

In the first step we read in the linelist. The one used in the example is that reported in the JQSRT 270 (2021) 107684. However this linelist is for demonstration purposes only, users should refer to that paper for the actual line list.

.. code:: ipython3

	from MATS.linelistdata import linelistdata
	PARAM_LINELIST = linelistdata['Singlet_Delta_Linelist_JQSRT_270_2021_107684']
	PARAM_LINELIST.sort_values('nu', inplace = True)
	
In this example, we have written a definition that simulates the monomer absorption for a spectrum, simulates a theoretical CIA based on the Karman et al model using the o2_cia_karman_model definition, and then adds the theoretical cia to the simulated monomer absorption in the spectrum alpha data.

The o2_cia_karman_model function requires wavenumber axis, temperature, pressure, and sample composition (O2 and N2 only), in addition to values for the Spin_orbit magnitude for both O2-O2 and O2-N2, exchange magnitude of O2-O2, temperature dependences for spin orbit (O2-O2 and O2-N2, can be constrained later to be equal) and exchange intensities. Additionally, the band being studied a_band or singlet_delta needs to be specified.


.. code:: ipython3

	wave_range = 1.5 #range outside of experimental x-range to simulate
	IntensityThreshold = 1e-30 #intensities must be above this value to be simulated
	Fit_Intensity = 1e-20 #intensities must be above this value for the line to be fit
	order_baseline_fit = 0


	wave_min = 7500 #cm-1
	wave_max = 8350 #cm-1
	wave_space = 0.01 #cm-1
	wavenumbers = np.arange(wave_min, wave_max, 0.02)


	def sim_spectra_with_CIA(pressure, temperature, sample_molefraction, filename):
		Diluent = {'O2': {'composition':sample_molefraction[7], 'm': 31.998}, 'N2': {'composition':1-sample_molefraction[7], 'm': 28.0134}}
		spec =  MATS.simulate_spectrum(PARAM_LINELIST, wavenumbers = wavenumbers, 
							 temperature = temperature,  pressure = pressure,
							 molefraction = sample_molefraction, Diluent = Diluent,
							 filename = filename, wing_wavenumbers = 25, wing_method = 'wing_wavenumbers',)
		# Adjusted from Initial guess derived based HITRAN 2020 reported theoretical CIA
		EXCH_c, EXCH_b, EXCH_a = [3.63e-06, 0.003, 1]
		SO_c, SO_b, SO_a =[1.5e-06, 0.0002, 1]
		SO_O2, SO_N2, EXCH_O2 = [40, 75 ,315] 

		CIA = MATS.o2_cia_karman_model(wavenumbers, temperature + 273.15, pressure/760, Diluent,
								SO_O2, SO_N2, EXCH_O2, 
								EXCH_b, EXCH_c, #O2-O2
								SO_b, SO_c, #O2-N2
								SO_b, SO_c, 
								SO_shift_O2_O2 = 0, SO_shift_O2_N2 = 0, EXCH_shift = 0,
								band = 'singlet_delta')
		spec.alpha += CIA
		return spec

	spec_1 = sim_spectra_with_CIA(760, 22.85, {7 :0.2095}, 'Air_296')
	spec_2 = sim_spectra_with_CIA(760, 0, {7 :0.2095}, 'Air_273')
	spec_3 = sim_spectra_with_CIA(760, 50, {7 :0.2095}, 'Air_323')
	spec_4 = sim_spectra_with_CIA(760, 22.85, {7 :1}, 'O2_296')
	spec_5 = sim_spectra_with_CIA(760, 0, {7 :1}, 'O2_273')
	spec_6 = sim_spectra_with_CIA(760, 50, {7 :1}, 'O2_323')


	spec_1.plot_wave_alpha()
	
.. image:: example_files/O2_CIA_simualted.png
							   
	
Set-up for Fitting
++++++++++++++++++
he CIA model is specified in the instantiation of the Dataset class. The CIA_model is specified through a dictionary specifying both the model and band. Currently, the Karman model is the only model available and the available bands are a_band and singlet_delta.

Like the generate_baseline_paramlist() definition, the generate_CIA_paramlist() generates a dataframe and file that summarizes the parameters based on the MATS preset values from the literature. https://doi.org/10.1038/s41557-018-0015-x, https://doi.org/10.1016/j.icarus.2019.02.034


.. code:: ipython3
	
	#Add all spectrum to a Dataset object
	SPECTRA = MATS.Dataset([spec_1, spec_2, spec_3, spec_4, spec_5, spec_6], 
						   'Experimental Singlet Delta with CIA', PARAM_LINELIST, CIA_model = {'model':'Karman', 'band': 'singlet_delta'})

	#Generate Baseline Parameter list based on number of etalons in spectra definitions and baseline order
	BASE_LINELIST = SPECTRA.generate_baseline_paramlist()
	CIA_LINELIST = SPECTRA.generate_CIA_paramlist()
	CIA_LINELIST
	
In the Generate_FitParam_File class, the CIA_linelist is set to the CIA_LINELIST generated by the generate_CIA_paramlist. The genererate_fit_KarmanCIA_linelist() will then add the uncertainty and boolean vary columns, setting the vary column to vary the parameters indicated in the definition call

.. code:: ipython3
	FITPARAMS = MATS.Generate_FitParam_File(SPECTRA, PARAM_LINELIST, BASE_LINELIST, CIA_linelist = CIA_LINELIST, 
											lineprofile = 'SDNGP', linemixing = True, 
									  fit_intensity = Fit_Intensity, threshold_intensity = IntensityThreshold, sim_window = wave_range,
									  nu_constrain = True, sw_constrain = True, gamma0_constrain = True, delta0_constrain = True, 
									   aw_constrain = True, as_constrain = True, 
									   nuVC_constrain = True, eta_constrain =True, linemixing_constrain = True)

	FITPARAMS.generate_fit_param_linelist_from_linelist(vary_nu = {7:{1:False, 2:False, 3:False}}, vary_sw = {7:{1:False, 2:False, 3:False}},
														vary_gamma0 = {7:{1: False, 2:False, 3: False}, 1:{1:False}}, vary_n_gamma0 = {7:{1:False}}, 
														vary_delta0 = {7:{1: False, 2:False, 3: False}, 1:{1:False}}, vary_n_delta0 = {7:{1:False}}, 
														vary_aw = {7:{1: False, 2:False, 3: False}, 1:{1:False}}, vary_n_gamma2 = {7:{1:False}}, 
														vary_as = {}, vary_n_delta2 = {7:{1:False}}, 
														vary_nuVC = {7:{1:False}}, vary_n_nuVC = {7:{1:False}},
														vary_eta = {}, vary_linemixing = {7:{1:False}})

	FITPARAMS.generate_fit_baseline_linelist(vary_baseline = False, vary_molefraction = {7:False, 1:False}, vary_xshift = False, 
										  vary_etalon_amp= False, vary_etalon_period= False, vary_etalon_phase= False, 
											 vary_pressure = False, vary_temperature = False)
	FITPARAMS.generate_fit_KarmanCIA_linelist(vary_S_SO = True, vary_S_EXCH = True, 
											vary_EXCH_temp = True, vary_SO_temp = True, 
											vary_EXCH_shift = False, vary_SO_shift = False)

Fit the SPECTRA
+++++++++++++++
Finally, in the Fit_Dataset instance, the CIA_linelist_file needs to be specified to indicate what parameters are being adjusted, provide initial values, and to allow uncertainties to be updated.

The CIA parameters can be constrained using the constrained_CIA parameters() definition. The default option is that the the spin orbit intensity temperature dependence is the same for O2-O2 and O2-N2. Additionally, the shift of the spin-orbit mechanism is normally constrained to be the same for O2-O2 and O2-N2. However, either of these constraints can be removed by setting S_temperature_dependence_constrained = False or shift_constrained = False.

The rest of the fitting runs in the same manner as for monomer-absorption only MATS fits.

.. code:: ipython3

	fit_data = MATS.Fit_DataSet(SPECTRA,'Baseline_LineList', 'Parameter_LineList', CIA_linelist_file = 'CIA_LineList',
								minimum_parameter_fit_intensity = Fit_Intensity)
	params = fit_data.generate_params()
	params = fit_data.constrained_CIA(params)
	for param in params:
		if 'S_SO' in param:
			params[param].set(min = 0.00)
		if 'S_EXCH' in param:
			params[param].set(min = 0.00)

	result = fit_data.fit_data(params, wing_wavenumbers = 25,  wing_method = 'wing_wavenumbers')
	fit_data.residual_analysis(result, indv_resid_plot=False)
	fit_data.update_params(result)
	SPECTRA.generate_summary_file(save_file = True)
	SPECTRA.plot_model_residuals()

.. image:: example_files/O2_CIA_fit.png