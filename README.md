# foxsi-response-model

Code for modeling the FOXSI instrument response. Descriptions of the various files are included below. 

### Main response modeling procedure files

The plots made in the following require some utility routines found in FOXSI SCIENCE, so to use them you should either have
that code in your IDL path (https://github.com/foxsi/foxsi-science), or you should re-write the plotting.

**response_modeling.pro** – Includes main response-modeling wrapper procedure, and a number of component functions used in the process.

**response_modeling_2thresh.pro** – Very similar to the above, but designed to work while using two thresholds (e.g. a trigger threshold 
and a lower "companion events" threshold used to determine whether there is a signal in the strips adjacent to the main trigger).  This is
the version which was used in the FOXSI ALS SPIE paper https://doi.org/10.1117/12.2629443, as well as in Jessie Duncan's dissertation 
https://hdl.handle.net/11299/241752.

**howditdo.pro** – How did it do? Simple, hardcoded method for comparing input sources and resulting sources (spatially) for the two examples 
used in the SPIE paper. Fits peaks and compares their locations, amplitudes, widths. Also makes a little comparison plot for visual reference.

### Other IDL procedures

**badpar.pro, lclxtrem.pro** – helpful utilities, Written by Marc W. Buie, Lowell Observatory. Included so you don't have to download the whole Buie library (here, if interested: https://www.boulder.swri.edu/~buie/idl/)

### Saved values for use in flare_spectrum function (to add simulated input M1 flare spectrum)
**m1_photon_spec.csv** - saved photon spectrum (M1 flare) from Yixian Zhang. 

**atten_factor.csv** - saved attenuation efficiency coeficients from Yixian Zhang. 

  Attenuation factors are used according to the following:
        We define, I = Io*EXP(-µx), the transmitted intensity (observed by detector) as a function of energy (X). 
        This assumes Io is the initial intensity (flare photon spectrum, e.g. values from m1_photon_spec.csv), and
        EXP(-µx) is the "attenuation efficiency" (values in atten_factor.csv, as a function of photon energy).
