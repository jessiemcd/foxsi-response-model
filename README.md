# foxsi-response-model

Code for modeling the FOXSI instrument response. Descriptions of the various files are included below. 



Saved values from Yixian Zhang, for simulated M1 flare spectrum. Read in by flare_spectrum function for use in response modeling.
m1_photon_spec.csv - saved photon spectrum (M1 flare)
atten_factor.csv - saved attenuation efficiency coeficients: 
                      I = Io*EXP(-µx) = transmitted intensity (observed by detector) as a function of energy (X)
                      Io = initial intensity (flare photon spectrum)
                      EXP(-µx) = "attenuation efficiency" (values in atten_factor.csv, as a function of photon energy)
