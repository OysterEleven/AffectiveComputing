.VoiceSource/

  * excitation_instants.m
	-- Finds the instances of excitation (glottal closure) in a
	   speech signal

  * closed_phases_adaptive.m
	-- Finds the closed phase region in a speech signal for each 
           glottal cycle

  * voiced_excite.m
	-- Computes a structure with voiced excitation information.

  * glottmarks.m
	-- Segments a speech signal into phonatory cycles and identifies
 	   glottal landmarks within each

  * LF.m
	-- Computes the glottal volume velocity waveform, and derivative,
	   using the Liljencrants-Fant parametric model, given the parameters.

  * fitLF.m
	-- Estimates the parameters of the LF model from glottal volume
	   velocity waveform data. Needs Optimization Toolbox.

  * glottal_invfilt.m
  	-- High-level function to identify phonatory cycles in a speech 
	   signal, inverse-filter and fit the parametric LF model to each 
           cycle.
  * gne.m
	-- Computes the glottal-to-noise excitation ratio from speech  

  * psp.m	
	-- Computes the parabolic spectral parameter from a cycle of the
	   glottal volume velocity waveform

  * jitt_shimm.m
	-- Computes jitter and shimmer values from speech.

.F0_HarmonicAnalysis/

  * getF0.m      
	-- Pitch extraction algorithm

  * dissmeasure.m
	-- Calculates the intrinsic dissonance value for a spectral pattern

  * timbre2disscurve.m   
	-- Calculates a dissonance curve for a spectral pattern

  * dissgram.m 
	-- Function to compute a perceptual dissonance diagram 
           (dissonance curves as a function of time) 

  * XtractHarmDissfeats.m 
	-- High-level function to extract features from a dissonance diagram. 


.Loudness/

  * loudnessgram.m
	-- Computes a loudness diagram (distribution of loudness value
	   over the Bark scale as a function of time).

  * XtractLoudnessfeats.m
	-- High-level function to extract features from a loudness diagram.

.Utils/

  * dapw.m
	-- Estimates parameters of a (frequency-weighted) discrete
	   all-pole model.

  * downsample_to_n_hz.m
	-- Wrapper for RESAMPLE() with extra error-handling.

  * find_flattest_seg.m
	-- Utility function to find flattest segment in a signal. Used to
	   help identify closed-phase region in glottal cycle by tracking
 	   first-formant movement.
  * find_peak.m
	-- General-purpose peak finder.

  * formant1a.m
	-- Estimates first formant from speech signal using the covariance
	   method.

  * formant_reshape.m
	-- "Fixes" an AR polynomial by discarding roots that meet
	   certain criteria. Type HELP FORMANT_RESHAPE for more details.

  * lpc_residual.m
	-- Estimates an AR model for a speech signal, inverse-filters, and
	   finds residual.

  * mypickpeak.m
	-- Wrapper for pickpeak.m to find peaks meeting certain criteria.

  * perturb.m
	-- Calculates the "perturbation factor" and "perturbation quotient"
	   of a sequence.
  * rms.m
	-- Finds the root-mean-square value of a sequence.

  * stft.m
	-- Special-purpose wrapper for a short-time Fourier transform 
	used to find excitation instants.

.Utils/DIN_2045631_V201.4.20/
		
	* DIN45631.m
		-- Calculate loudness based on DIN 45631 standard.

	* Oct3dsgn.m
		-- Designs bank of 1/3-octave filters.

	* myGenerateFilters.m
		-- Wrapper to generate 1/3-octave filters with OCT3DSGN().

.Utils/Voicebox

	* lpcar2fm.m
		-- Converts AR coefficients to formant freq+amp+bw triplets.

	* lpcar2zz.m
		-- Converts AR coefficients to z-plane poles.

	* lpcauto.m
		-- Estimates AR coefficients using autocorrelation analysis.

	* lpccovar.m
		-- Estimates AR coefficients using covariance analysis.

	* lpcifit.m
		-- Inverse filters a speech signal given an AR model.

	* window.m
		-- Generates a standard windowing function.

	* windowx.m
		-- Generates ordinates for WINDOW().	
