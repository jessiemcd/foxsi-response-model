FUNCTION flare_spectrum

	;of form [energy (keV), mass attenuation coefficient of Al (cm^2/g)]
	;units cm^2/g
	;From: https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z13.html
	mass_atten = [[4., 3.605E+02], [5., 1.934E+02], [6., 1.153E+02], [8., 5.033E+01], $
					[10., 2.623E+01], [15., 7.955E+00], [20., 3.441E+00], [30., 1.128E+00], $
					[40., 5.685E-01]]

	coeffs = mass_atten
	for i=0, n_elements(mass_atten)/2-1 do begin
		;density of aluminum: 2.7 g/cm^3
		;using 260 um aluminum (0.026 cm)
		pair = mass_atten[*,i]
		replen = 2.7*pair[1]*0.01
		coeff = EXP(-replen)
		coeffs[1,i] = coeff
	endfor

	print, coeffs
	;plot, coeffs[0,*], coeffs[1,*]

	spectrum = read_csv('m1_photon_spec.csv') 
	atten = read_csv('atten_factor.csv') 
	
	atten_interp = INTERPOL(atten.field2, atten.field1, spectrum.field1)
	
	atten_spectrum = spectrum.field2*atten_interp
	
	;plot, atten.field1, atten.field2, xrange=[0,45], yrange=[0,0.1]
	;oplot, spectrum.field1, atten_interp, linestyle=2
	
	;plot, spectrum.field1, spectrum.field2, /ylog
	;oplot, spectrum.field1, atten_spectrum, linestyle=2

	

	;spectrum.field1
	;spectrum.field2

	rf2 = round(atten_spectrum)

	samples=[]
	for i=300, n_elements(spectrum.field1)-1 do begin

		if rf2[i] GT 0 then begin
			;This is for using the more sparse attenuation array from NIST mass attenuation coeffs
;			if spectrum.field1[i] GT 4 and spectrum.field1[i] LE 5 then num = rf2[i]*coeffs[1,0]
;			if spectrum.field1[i] GT 5 and spectrum.field1[i] LE 6 then num = rf2[i]*coeffs[1,1]
;			if spectrum.field1[i] GT 6 and spectrum.field1[i] LE 8 then num = rf2[i]*coeffs[1,2]
;			if spectrum.field1[i] GT 8 and spectrum.field1[i] LE 10 then num = rf2[i]*coeffs[1,3]
;			if spectrum.field1[i] GT 10 and spectrum.field1[i] LE 15 then num = rf2[i]*coeffs[1,4]
;			if spectrum.field1[i] GT 15 and spectrum.field1[i] LE 20 then num = rf2[i]*coeffs[1,5]
;			if spectrum.field1[i] GT 20 and spectrum.field1[i] LE 30 then num = rf2[i]*coeffs[1,6]
;			if spectrum.field1[i] GT 30 and spectrum.field1[i] LE 40 then num = rf2[i]*coeffs[1,7]
;			if spectrum.field1[i] GT 40 then num = rf2[i]
			num=rf2[i]
			arr = make_array(num, val=spectrum.field1[i])
			samples = [samples, arr]
		endif
		;print, i

	endfor
	print, size(samples)
	
	write_csv, 'm1_photon_list.csv', samples
	;plot_hist, samples, bin=0.01, /oplot, linestyle=3
	
	
	return, samples

END

;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================

FUNCTION add_sharing_new, distribution=distribution, threshold=threshold, $
			input_spectrum=input_spectrum, method=method, erez=erez, edist=edist, $
			stripbins=stripbins, side=side, sec_threshold=sec_threshold

	ps = size(distribution)
	
	number=5000.
	if edist EQ 'ALS' then begin
		es1 = randomu(seed, number*25./31., /normal) + 7
		es2 = randomu(seed, number*5./31., /normal) + 21
		es3 = randomu(seed, number/31., /normal) + 28
		es = [es1, es2, es3]
		se = size(es)
	endif
	
	if edist EQ 'crudeAm241' then begin
		;It's called crude Am241 bc only a vague effort has been made to get the line
		;ratios anywhere close to reality
		es1 = randomu(seed, number*10./74., /normal) + 13.9
		es2 = randomu(seed, number*15./74., /normal) + 18
		es3 = randomu(seed, number*5/74., /normal) + 21
		es4 = randomu(seed, number*4/74., /normal) + 26
		es5 = randomu(seed, number*40/74., /normal) + 60
		es = [es1, es2, es3, es4, es5]
		se = size(es)
	endif
	
	if edist EQ 'M1' then begin
		flare = read_csv('m1_photon_list.csv') 
		samples = flare.field1
		ss = size(samples)
		
		;plot_hist, samples, bin=0.5
	endif
	
	if edist EQ 'ALS_high' then begin
		es = randomu(seed, number, /normal) + 21
		se = size(es)
	endif
	


	;Taking the convolved source + converting to an event list with energies
	;For each bin, round to the nearest integer number of counts
	;For each count, add a new event to the list like: [bin #, energy]
	events_list=[]
	input_spectrum=[]
	ee_s=[]
	inds=[]
	for c=0, n_elements(distribution)-1 do begin
		for e=0, round(distribution[c]) do begin
			;Energies randomly generated between 0-100 ("keV")
			if edist EQ 'uniform' then ee = randomu(seed)*100
			;Energies set to equal pre-made distribution (see above)
			if edist EQ 'ALS' or edist EQ 'crudeAm241' or edist EQ 'ALS_high' then begin
				ind = fix(randomu(seed)*(se[1]-1))
				ee = es[ind]
				;print, ee
			endif
			;Energies set to equal pre-made imported distribution (M1 flare spectrum)
			if edist EQ 'M1' then begin
				ind = round(randomu(seed, /double)*(ss[1]))
				while ind GT ss[1]-1 do begin
					ind = round(randomu(seed, /double)*(ss[1]))
				endwhile
				ee = samples[ind]
				ee_s = [ee_s, ee]
				inds = [inds, ind]
			endif
			
			;ee = 50.
			input_spectrum = [input_spectrum, ee]
			events_list = [[events_list], [c, ee]]
		endfor
	endfor

	ss = size(events_list)
	print, 'Size of Events List: ', ss
	
	
	if method EQ 'alssharing' then begin
		;size of central region: (must be an integer # of micron, as we are implementing 1um bins)
		if side EQ 'Pt' then begin
			censize=round(44)
		endif
		
		if side EQ 'Al' then begin
			censize=round(20)
		endif
	
		;NOT DONE HERE- WORKING TO DEFINE THREE ZONES, NEED TO MAKE SURE THEY ALL COME OUT TO BE INTEGER BINS
		cross_wid = (stripbins - censize)/2
		M = cross_wid
	endif
		
	

	counter=0
	frame_list = []
	ogbins = []
	energy_ratios = []
	E1slost=[]
	E2slost=[]
	lostcomponents=[]
	keptcomponents=[]
	eventslost=[]
	eventslost1=[]
	;For every event in the event list:
	for ei=0, ss[2]-1 do begin
		;Get pair (location, energy)
		entry = events_list[*,ei]
		ogbins = [ogbins, (entry[0] mod fix(stripbins))]
		;Find the strip that contains the location
		hit_strip = fix(entry[0])/fix(stripbins)
		frame = make_array(2,2,val=0.)
		
		;If the location is greater than halfway through the first strip and less than one 
		;strip-half from the end of the array then continue (we don't want sharing with "strips" 
		;outside our dataset/detector)
		if entry[0] GE stripbins/2. and entry[0] LT ps[1]-stripbins/2. then begin
				
				if method EQ 'alssharing' then begin 
					;This method was developed later (spring 2022), and does not have the same saved
					;components, etc. for checking performance as with the other methods. 
					
						;If the event position is in the central region, assign all energy to the hit strip.
						if ogbins[-1] GE cross_wid and ogbins[-1] LT cross_wid+censize then begin
		
							;FOXSI ENERGY RESOLUTION ESTIMATE AT 1.2 KEV
							;SIGMA = 1.2/2.335 -> SIGMA=~0.51
			
							if erez EQ 'erez_on' then begin
								change = randomu(seed, /normal)*0.51
								if entry[1]+change GT threshold then frame[*,0] = [hit_strip, entry[1]+change]
							endif
			
							if erez EQ 'erez_off' then begin
								if entry[1] GT threshold then frame[*,0] = [hit_strip, entry[1]]
							endif
							;Add frame to list of frames and continue to the next event/frame
							frame_list = [[[frame_list]], [[frame]]]
						endif
		
						;If the event position is in the charge sharing region, do that.
						if ogbins[-1] LT cross_wid or ogbins[-1] GE cross_wid+censize then begin
		
							;If the location is greater than bin 5 and more than 5 bins from the end of the 
							;array then continue (we don't want sharing with "strips" outside our dataset/detector)
							if entry[0] GE cross_wid and entry[0] LT ps[1]-cross_wid then begin
									;Depending on which side of the strip the location falls, define the adjacent
									;strip (either the one above or below)
									next_strip = round((fix(entry[0]) mod fix(stripbins))/float(stripbins))
									if next_strip EQ 0 then begin
										next_strip = -1
										boundind = hit_strip*fix(stripbins)
									endif
									if next_strip EQ 1 then boundind = (hit_strip+1)*fix(stripbins)
									next_strip+=hit_strip
					
									;NEXT: ADD ENERGY ASSIGNMENT MATH!!
									S1 = min([hit_strip, next_strip])
									S2 = max([hit_strip, next_strip])
					
									E1 = entry[1]/2 * (1+(boundind - entry[0])/M)
									E2 = entry[1]/2 * (1-(boundind - entry[0])/M)

									;SANITY CHECK
									if E1 LT 0 or E2 lt 0 then STOP
			
									;Add strip numbers to frame
									frame[0,0] = S1
									frame[0,1] = S2
				
									;Add strip #, energy pairs (with EREZ) to each frame ONLY IF ENERGY ABOVE THRESHOLD
									if erez EQ 'erez_on' then begin
										change = randomu(seed, /normal)*0.51
										nE1 = E1+change 
										change = randomu(seed, /normal)*0.51
										nE2 = E2+change 
										
										if nE1 GT threshold then begin
											frame[1,0] = nE1
											if nE2 GT sec_threshold and nE2 LE threshold then frame[1,1] = nE2
										endif
										if nE2 GT threshold then begin
											frame[1,1] = nE2
											if nE1 GT sec_threshold and nE1 LE threshold then frame[1,0] = nE1
										endif
									endif
					
									;Add strip #, energy pairs (without EREZ) to each frame ONLY IF ENERGY ABOVE THRESHOLD
									if erez EQ 'erez_off' then begin
										if E1 GT threshold then begin
											frame[1,0] = E1
											if E2 GT sec_threshold and E2 LE threshold then frame[1,1] = E2
										endif
										if E2 GT threshold then begin
											frame[1,1] = E2
											if E1 GT sec_threshold and E1 LE threshold then frame[1,0] = E1
										endif
									endif
									;Add frame to list of frames and continue to the next event/frame
									frame_list = [[[frame_list]], [[frame]]]
							endif
						endif
				endif
		endif
		if 	entry[0] LT stripbins/2. or entry[0] GE ps[1]-stripbins/2. then begin
			hit_strip = fix(entry[0])/fix(stripbins)
			frame = make_array(2,2,val=0.)
			if entry[1] GT threshold then frame[*,0] = [hit_strip, entry[1]]
			;frame_list = [[[frame_list]], [[frame]]]	
		endif
	endfor
	

	;save, frame_list, filename='fakedata.sav'
	print, 'Frame List Size:', size(frame_list)
	;print, 'Input Spectrum Size:',size(input_spectrum)
	
	
	save, input_spectrum, filename='input_spectrum.sav'
	
	return, frame_list
END




;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================

FUNCTION process_frame_list, frame_list=frame_list
	s = size(frame_list)

	raw_spectrum=[]
	strip_list=[]
	;For each frame in the frame list
	for i=0, s[3]-1 do begin
		;Extract the specific frame (2x2)
		f = frame_list[*,*,i]
		;Energies from the frame
		es = f[1,*]
		;Find where the energies are non-zero + add non-zero energies to spectrum
		register = where(es GT 0)
		if register[0] NE -1 then raw_spectrum = [raw_spectrum, es[register]]
		;Add strips with non-zero energies to list of strips
		ss = f[0,*]
		if register[0] NE -1 then strip_list=[strip_list, ss[register]]
	endfor
	print, 'Raw Spectrum Size: ', size(raw_spectrum)
	print, 'Strip List Size: ', size(strip_list)
	return, [[[raw_spectrum]], [[strip_list]]]
END

;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================

FUNCTION remove_sharing_new, frame_list=frame_list, method=method, erez=erez, $
			threshold=threshold, sec_threshold=sec_threshold, $
			stripbins=stripbins, side=side, siz=siz
	;What do we want to do here?
	
	;Take the input frame list (1-2 registered strips with energies)
	;Revert back to 1 event per frame:
		;- For single-strip events, randomly assign position within strip
		;- For double-strip events, use (modified) expression from furukawa (2020) to
		;	extract a position based on the two strips and two energies registered. 
		
		;We should be back to the same number of above-threshold events we had originally 
		; (we lose events where the total original energy is less than the threshold used in 
		;	add_sharing). 

	s = size(frame_list)
	
	
	if method EQ 'alssharing' then begin
	
		if stripbins NE 60 then begin
			print, 'ALS STYLE CHARGE SHARING ONLY IMPLEMENTED FOR 60 BINS PER STRIP'
			PRINT, 'YOU ASKED TO DO IT WITH ', stripbins
			stop
		endif
				
		if side EQ 'Pt' then begin
			yellow = 5./50.
			censize=round(44)
		endif
		
		if side EQ 'Al' then begin
			yellow = 12./32.
			censize=round(20)
		endif
	
		;size of central region: (must be an integer # of micron, as we are implementing 1um bins)
		cross_wid = (stripbins - censize)/2
		M = cross_wid
	endif
		
	
	spectrum=[]
	hit_list=[]
	strip_list=[]
	double_bins = []
	single_bins=[]
	rawcsa_spectrum =[]
	rec_comp_spectrum = []

	
	addedback=[]
	bin_value1=[]
	bin_value2=[]
	
	bin_dub=[]
	
	;For each frame in the frame list
	for i=0, s[3]-1 do begin
		;Extract the specific frame (2x2)
		f = frame_list[*,*,i]
		;Energies from the frame
		es = f[1,*]
		
		
		register = where(es GT 0)
		
		;Add strips with non-zero energies to list of strips
		ss = f[0,*]
		if register[0] NE -1 then strip_list=[strip_list, ss[register]]
		if register[0] NE -1 then rawcsa_spectrum = [rawcsa_spectrum, total(es)]
		
		;If there's only a single event in the frame:
		if n_elements(register) EQ 1 and register[0] NE -1 then begin
			
			if method EQ 'alssharing' then begin
			
				;If a random number between 0-1 is greater than the given number (% of events with loss)
				;add back a below-threshold component in one of the two adjacent strips 
				
				;Setting the following condition to be impossible (I guess EXTREMELY unlikely) as not using 
				;the add-back with the two-threshold method - although did slightly widen single-strip 
				;assignment region, see below. response
				if randomu(seed, 1) LE -100 then begin
					print, 'Doing this'

					if f[1,1] EQ 0 then begin
						;while f[1,1] LE 0. or f[1,1] GT 4. do begin
						while f[1,1] LE 0. or f[1,1] GT sec_threshold do begin
							newe = randomu(seed)*sec_threshold
							;newe = randomu(seed)*4.
							f[1,1] = newe
							addedback = [addedback, newe]
						endwhile
						which = RANDOMU(SEED, /NORMAL) GE 0
						if which EQ 0 then f[0,1] = f[0,0]-1
						if which EQ 1 then f[0,1] = f[0,0]+1
					endif

					if f[1,0] EQ 0 then begin
						;while f[1,0] LE 0. or f[1,0] GT 4. do begin
						while f[1,0] LE 0. or f[1,0] GT sec_threshold do begin
							newe = randomu(seed)*sec_threshold
							;newe = randomu(seed)*4.
							f[1,0] = newe
							addedback = [addedback, newe]
						endwhile
						which = RANDOMU(SEED, /NORMAL) GE 0
						if which EQ 0 then f[0,0] = f[0,1]-1
						if which EQ 1 then f[0,0] = f[0,1]+1
					endif 

					;switch which event is where if their indices are out of order
					if f[0,0] GT f[0,1] then begin
						newframe = make_array(2,2,val=0.)
						newframe[*,0] = f[*,1]
						newframe[*,1] = f[*,0]
						f = newframe
					endif

					es = f[1,*]
					;re-define register so this gets caught as a double strip event later.
					register = where(es GT 0)
	;;				
				endif else begin

					;For the single strip events we are actually going to treat like single strip
					if f[1,1] EQ 0 then strip = f[0,0]
					if f[1,0] EQ 0 then strip = f[0,1]
				
					;Note: to reduce periodicity from "gaps" due to events with a below-threshold component,
					;center region where single-strip events are placed has been slightly widened. 
					bin = cross_wid-2 + round((censize+4)*randomu(seed, 1))
					;bin = cross_wid + round((censize)*randomu(seed, 1))
					if bin GT 60 then print, 'Huge Bin:', bin
					X = strip*60 + bin
					bin_value1 = [bin_value1, bin]
					if X LT siz and X GE 0 then begin
						;Add reconstructed position to hit list + reconstructed energy to energy list
						hit_list = [hit_list, X]
						spectrum = [spectrum, total(f[1,*])]
					endif
				
				endelse
			
			endif
			
		endif
		
		;If there are two events in the frame
		if n_elements(register) EQ 2 then begin
			
			if method EQ 'alssharing' then begin
			
			
				Re = (f[1,0] - f[1,1])/total(f[1,*])
				boundind = max(f[0,*])*stripbins
				X = boundind - Re*cross_wid
				
				bin_dub = [bin_dub, Re*cross_wid]
			
				;print, 'Subtracted factor:', Re*cross_wid
			
				if X LT siz and X GE 0 then begin
					;Add reconstructed position to hit list + reconstructed energy to energy list
					hit_list = [hit_list, X]
					spectrum = [spectrum, total(f[1,*])]
				endif
			
				xbin = X mod fix(stripbins)
			
				if xbin gt stripbins then print,  'Huge bin:', xbin
				;if xbin lt 8 then print, xbin, f
				bin_value2 = [bin_value2, xbin]
			
			endif
			
		endif
		
		es = f[1,*]
		register = where(es GT 0)
		rec_comp_spectrum = [rec_comp_spectrum, es[register]]
		
	endfor
	
	;plot_hist, bin_dub, bin=1
	;print, min(bin_dub)
	;print, max(bin_dub)
	
	;stop

	
	save, spectrum, rawcsa_spectrum, rec_comp_spectrum, filename='saved_restored_spectrum.sav'
	
	rc_hist = histogram(hit_list, binsize=1, nbins=siz, min=0)
	return, rc_hist
	
END



;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================

function PROB, y, xi, stripbins, sigma=sigma, PSF_FWHM=PSF_FWHM
;(USING GAUSSIAN) - this is our optics PSF analogue.
;For a single y value: probability that the measured source location will fall
;into the interval (y, y+dy) when xi is the source location.
;For an array of y values: returns array of same size (probability of each y).
;Values in said array will sum to 1 if the source is located sufficiently far
;from edges of y array, to less if not.

  if KEYWORD_SET(PSF_FWHM) then sigma = PSF_FWHM/2.355/10*stripbins

  probb = (SQRT(2.*!PI)*sigma)^(-1)*EXP(-((y-xi)^2/(2.*sigma^2)))
  ;probb = make_array(n_elements(y), val=1./n_elements(y))
  return, probb
end



function fancy_PROB, y, xi, stripbins
;This is a better estimate of the optics PSF (still axially symmetric) based on parameters
;sent by Milo to describe a fit to the actual measured PSF of a FOXSI optic (X7, a 10-shell
;optic). 
;For a single y value: probability that the measured source location will fall
;into the interval (y, y+dy) when xi is the source location.
;For an array of y values: returns array of same size (probability of each y).
;Values in said array will sum to 1 if the source is located sufficiently far
;from edges of y array, to less if not.


A1 = 0.151
A2 = 0.383
A3 = 0.466
S1 = 6.677/10*stripbins
S2 = 1.111/10*stripbins
S3 = 2.580/10*stripbins

probb = A1*(SQRT(2.*!PI)*S1)^(-1)*EXP(-((y-xi)^2/(2.*S1^2))) + $
		A2*(SQRT(2.*!PI)*S2)^(-1)*EXP(-((y-xi)^2/(2.*S2^2))) + $
		A3*(SQRT(2.*!PI)*S3)^(-1)*EXP(-((y-xi)^2/(2.*S3^2)))
  ;probb = make_array(n_elements(y), val=1./n_elements(y))
  return, probb
end

;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================


FUNCTION ogpsf, xis, xs, stripbins, PSF_FWHM=PSF_FWHM
;ogpsf does the same thing as modpsf, except it does not include the
;pixelization function, pixelize.

  y = indgen(n_elements(xis))
  sigma=3.
  npsf=[]
  for i=0, n_elements(xis)-1 do begin
    col = []
    for j=0, n_elements(xs)-1 do begin
      ;val = PROB(xs[j], xis[i], sigma)
      val = fancy_PROB(xs[j], xis[i], stripbins)
      if keyword_set(PSF_FWHM) then val = PROB(xs[j], xis[i], stripbins, PSF_FWHM=PSF_FWHM)
      col=[col,val]
    endfor
    npsf = [[npsf], [col]]
  endfor
  S = size(npsf)
  IF S[1] GT 1 then begin
      ;If the psf is 2D, let's have a look at it!
      ct = colortable(1, /reverse, stretch=-60)
      ;g = image(npsf, AXIS_STYLE=2, MARGIN=0.1, xtitle='X Values', ytitle='XI Values', $
      ;    RGB_TABLE=ct, title='OG PSF')
  ENDIF
  return, npsf
END

;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================

FUNCTION deconvolve_simple, data=data, iter=iter, og=og, matrix=matrix
	;Currently hardcoded for data, source to have the same dimensions

	xis = findgen(n_elements(data))
	xs = findgen(n_elements(data))
	guess = make_array(n_elements(data), value=1.)

	psf_totals_xi = []
	for i=0, n_elements(xis)-1 do begin
		tot = total(matrix(*, xis[i]))
		psf_totals_xi = [psf_totals_xi, tot]
	endfor	
	;print, psf_totals_xi
	
	;setting initial psi^r to initial guess, then iteratively updating it
	deconvs=[]
	coors=[]
	psir_ = guess
	ii=0
	phirxs=[]
	while (ii LT iter) do begin

		phir_x=[]
		for i=0, n_elements(xs)-1 do begin
			ff = psir_*matrix(xs[i], *)
			phi = total(ff)
			phir_x=[phir_x, phi]
		endfor
	
		phirxs = [[phirxs], [phir_x]]

		psir1_xi=[]
		for i=0, n_elements(xis)-1 do begin
			rat = (data/phir_x)
			rat[where(finite(rat) EQ 0)] = 0.0
			;plot, rat
			ff = rat*matrix(*, xis[i])
			;intt = int_tabulated(xs, ff)
			intt = total(ff)
			psir_val = psir_[xis[i]]
			psii = psir_val*intt
		
			psir1_xi=[psir1_xi, psii]
		endfor
		;print, size(psir1_xi)

		;This step was introduced to reduced artificial periodicity when using a pixelated
		;method. It doesn't seem to do much when we're only using the gaussian PSF (all 1s 
		;except for at the far edges).
		;psir1_xi = psir1_xi/psf_totals_xi
	
		psir_ = psir1_xi
		ii+=1
		;print, size(psir1_xi)
		;print, size(og)
		coor = c_correlate(psir1_xi, og, 0)
		;print, 'Coor:', coor
		;print, og
		;print, psir1_xi
    	;coors = [coors, coor]
    	psir1_xi = [psir1_xi, coor]
		deconvs= [[deconvs], [psir1_xi]]
	endwhile

	return, deconvs
END




;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================
;========================================================================================


pro response_modeling_2thresh, source=source, strips=strips, threshold=threshold, $
						sec_threshold=sec_threshold, method=method, erez=erez, edist=edist, $
						stripbins=stripbins, edge=edge, iter=iter, $
						FULLPLOT=FULLPLOT, COMPACTPLOT=COMPACTPLOT, PLOTSPECTRA=PLOTSPECTRA, $
						PSF_FWHM=PSF_FWHM
						
						
;A copy of response_modeling.pro, in which we will edit the various relevant functions such that
;the applied threshold (4 keV) applies to the higher-energy component, with a lower 1.4 keV secondary
;threshold for adjacent strips. This is intended to be the method used in the SPIE proceedings paper 
;(Duncan et al. 2022, July SPIE Astronomical Telescopes + Instrumentation)

;Example syntax
;response_modeling_2thresh, source='twosource', method='alssharing', erez='erez_on', edist='ALS_high', edge=1, iter=60, /fullplot

;Currently hardcoded using Al-side parameters if ALS-style charge sharing (energy ratios) is 
;implemented.		
side = 'Al'

;KEYWORDS
;
; SOURCE - SOURCE TYPE (spatial) (options: 
;			twosource: two gaussians
;			twosource_close: same two gaussians, one strip width apart
;			uniform: uniform
;			faint10xbright: one gaussian 10x brighter than a second
;			bfsb: broad faint + small bright gaussians
;			pinslit: approzimating "pinhole"/slit shape from furukawa et al. 2020
;	
; EDIST - SOURCE TYPE (spectral) (options:
;			ALS: simulate 7 keV beam with 21, 28 keV harmonics
;			ALS_high: 21 keV beam only
;			crudeAm241: It's called crude Am241 bc only a vague effort has been made to get the line
;						ratios anywhere close to reality
;			M1: M1 flare spectrum simulated by Yixian (see flare_spectrum function for how 
;						this is prepared)
;			uniform: uniform spectral distribution (random energies between 0-100 keV)
;
; EDGE - are sources centered on the edges of strips or in their centers? (0=edge centers)
;
; STRIPS - how many strips to simulate
; STRIPBINS - how many bins per strip?
; THRESHOLD - detector threshold (energy space)
; METHOD - which charge sharing method to use. 
;				alssharing: Use ALS energy ratio sharing method (newer, circa spring 2022) 
;							ONLY METHOD INCLUDED HERE - response_modeling.pro has some others
;
; ITER - how many deconvolution iterations to do before stopping. 
; PSF_FWHM - if set, this keyword allows user-inputed PSF FWHM (gaussian PSF) in ~arcsec. The default is
;				to use the measured PSF from milo (FANCY-PROB)
;
;Plotting keywords - can set any or all of these. You should set at least one or this doesn't 
;                    show you results in a very clear way!
;
; FULLPLOT - set to make a six-panel plot with the original source, convolved source, strip image,
;				reconstructed image (pre-deconv), and deconvolution plot
;
; COMPACTPLOT - set to make a smaller plot with just the original source, strip image, and the
;				deconvolution plot
;
; PLOTSPECTRA - set to plot spectra



;What questions do we have and what distributions do we want to compare over the course
;of this process?

;Energy - Compare:
;generated fake spectrum 
;raw (charge shared) spectrum 
;reconstructed spectrum

;Position Information - Compare
;Initial distribution
;Convolved distribution
;Strip-binned (charge shared) distribution
;Reconstructed distribution
;Deconvolved reconstructed distribution


;add_sharing takes in: 
;-convolved distribution
;-threshold
;add_sharing returns
;-frame list
;
;remove_sharing takes in:
;-frame list
;remove_sharing returns:
;-raw spectrum
;-CSA spectrum
;-reconstructed distribution

default, edge, 0
default, strips, 20
default, stripbins, 60
default, iter, 20
default, source, 'twosource'
default, threshold, 4
default, sec_threshold, 1.4
default, method, 'alssharing'
default, erez, 'erez_on'
default, edist, 'uniform'

siz = strips*stripbins

;FOR PLOTTING IN POSITION SPACE: INDICES TO PLOT (E.G. FOR SUBSET OF TOTAL ARRAY)
indices=[5*stripbins:12*stripbins]
indices_strip = [5:12]


;Here we have defined a source and recorded data of the same size. 
xs = findgen(siz)
xis = findgen(siz)

;========================================================================================
;========================================================================================


if source EQ 'two_als_beams' then begin
	;Generating observed source: two 2 µm beams 6 µm apart - ONLY TO BE USED WITH 60 BINS/STRIP
	  loc1 = 479
	  loc2 = 482
	  sigma1 = 0.8
	  sigma2 = 0.8
	  mag1 = 100.
	  mag2 = 100.

	  vall=10.
	  psi_xi = mag1*EXP(-((xis-loc1)^2/(2.*sigma1^2))) + $
			  mag2*EXP(-((xis-loc2)^2/(2.*sigma2^2)))
	  psi_xi = psi_xi + make_array(n_elements(xis), val=vall)  
	  
	  print, source, ':'
	  print, 'Center 1: ', loc1, ' Center 2:', loc2
	  print, 'Width 1: ', sigma1, ' Width 2:', sigma2
	  save, loc1, loc2, sigma1, sigma2,  filename=source+'_inputs.sav'
endif 


if source EQ 'twosource' then begin
	;Generating observed source: two peaks on continuum
	  loc1 = 7.*stripbins
	  loc2 = 10.*stripbins
	  if edge EQ 0 then BEGIN
	  		loc1+=0.5*stripbins
	  		loc2+=0.5**stripbins
	  endif
	  sigma1 = 1./10*stripbins
	  sigma2 = 2./10*stripbins
	  mag1 = 50.
	  mag2 = 70.

	  vall=15.
	  psi_xi = mag1*EXP(-((xis-loc1)^2/(2.*sigma1^2))) + $
			  mag2*EXP(-((xis-loc2)^2/(2.*sigma2^2)))
	  psi_xi = psi_xi + make_array(n_elements(xis), val=vall)  
	  
	  print, source, ':'
	  print, 'Center 1: ', loc1, ' Center 2:', loc2
	  print, 'Width 1: ', sigma1, ' Width 2:', sigma2
	  save, loc1, loc2, sigma1, sigma2,  filename=source+'_inputs.sav'
endif 

if source EQ 'twosource_close' then begin
	;Generating observed source: two peaks on continuum close together
	
	  loc1 = 7.*stripbins
	  loc2 = 8.*stripbins
	  if edge EQ 0 then BEGIN
	  		loc1+=0.5*stripbins
	  		loc2+=0.5*stripbins
	  endif
	  sigma1 = 1./10*stripbins
	  sigma2 = 2./10*stripbins
	  mag1 = 50.
	  mag2 = 70.
	  
	  
	  vall=15.
	  psi_xi = mag1*EXP(-((xis-loc1)^2/(2.*sigma1^2))) + $
			  mag2*EXP(-((xis-loc2)^2/(2.*sigma2^2)))
	  psi_xi = psi_xi + make_array(n_elements(xis), val=vall)  
	  print, source, ':'
	  print, 'Center 1: ', loc1, ' Center 2:', loc2
	  print, 'Width 1: ', sigma1, ' Width 2:', sigma2
	  save, loc1, loc2, sigma1, sigma2,  filename=source+'_inputs.sav'
endif 

if source EQ 'faint10xbright' then begin
	;Generating observed source: two peaks on continuum, one 10x brighter than the other
	  loc1 = 7.*stripbins
	  loc2 = 10.*stripbins
	  if edge EQ 0 then BEGIN
	  		loc1+=0.5*stripbins
	  		loc2+=0.5*stripbins
	  endif
	  sigma1 = 1./10*stripbins
	  sigma2 = 2./10*stripbins
	  mag1 = 60.
	  mag2 = 600.
	  
	  vall=10.
	  psi_xi = mag1*EXP(-((xis-loc1)^2/(2.*sigma1^2))) + $
			  mag2*EXP(-((xis-loc2)^2/(2.*sigma2^2)))
	  psi_xi = psi_xi + make_array(n_elements(xis), val=vall)  
	  print, source, ':'
	  print, 'Center 1: ', loc1, ' Center 2:', loc2
	  print, 'Width 1: ', sigma1, ' Width 2:', sigma2
	  save, loc1, loc2, sigma1, sigma2,  filename=source+'_inputs.sav'
endif 

if source EQ 'bfsb' then begin
	;Generating observed source: small bright source on top of broad faint source
	  loc1 = 7.*stripbins
	  loc2 = 9.*stripbins
	  if edge EQ 0 then BEGIN
	  		loc1+=0.5*stripbins
	  		loc2+=0.5*stripbins
	  endif
	  sigma1 = 1./10*stripbins
	  sigma2 = 15./10*stripbins
	  mag1=40.
	  mag2=20.

	  
	  psi_xi = mag1*EXP(-((xis-loc1)^2/(2.*sigma1^2))) + $
			  mag2*EXP(-((xis-loc2)^2/(2.*sigma2^2)))
			  
		print, source, ':'
	  print, 'Center 1: ', loc1, ' Center 2:', loc2
	  print, 'Width 1: ', sigma1, ' Width 2:', sigma2
	  save, loc1, loc2, sigma1, sigma2,  filename=source+'_inputs.sav'
endif 


if source EQ 'pinslit' then begin
	;Generating observed source: single peak with ~100µm width
	  loc2 = 9.*stripbins
	  if edge EQ 0 then BEGIN
	  		loc2+=0.5*stripbins
	  endif
	  sigma2 = 8./10*stripbins
	  mag2=50.
	  
	  psi_xi = mag2*EXP(-((xis-loc2)^2/(2.*sigma2^2)))
endif 


if source EQ 'uniform' then begin
	;Generating observed source: uniform counts/bin
	  psi_xi = make_array(n_elements(xis), val=100.) 
endif 




;========================================================================================
;========================================================================================

og_source = psi_xi

;plot, og_source

print, 'OG total counts:', total(og_source)

;CONVOLUTION - simply the product of the PSF matrix and the source array
matrix = ogpsf(xis, xs, stripbins, PSF_FWHM=PSF_FWHM)

phi_psf = matrix_multiply(matrix, psi_xi)
;oplot, phi_psf

print, 'Convolved total counts:', total(phi_psf)

input_spectrum=[]

frame_list = add_sharing_new(distribution=phi_psf, threshold=threshold, sec_threshold=sec_threshold, $
				input_spectrum=input_spectrum, erez=erez, $
				edist=edist, side=side, stripbins=stripbins, method=method)


print, 'Frames: ', size(frame_list)


pfl = process_frame_list(frame_list=frame_list)

raw_spectrum = pfl[*,*,0]
strip_list = pfl[*,*,1]

sl_hist = histogram(strip_list, binsize=1, nbins=20, min=0)
;plot, sl_hist, xtitle='X (Strips)', ytitle = 'Counts', thick=th, $
;		title='Strip-binned Source (Post Charge Sharing)',linestyle=6, psym=4

rcs = remove_sharing_new(frame_list=frame_list, method=method, erez=erez, $
						threshold=threshold, sec_threshold=sec_threshold, $
						stripbins=stripbins, side=side, siz=siz)
;plot, rcs

;averaging results - the idea is that this smooths out noise. Arbitrarily averaging 10 times
averages = 10
while averages GT 0 do begin
	rcs2 = remove_sharing_new(frame_list=frame_list, method=method, erez=erez, $
						threshold=threshold, sec_threshold=sec_threshold, $
						stripbins=stripbins, side=side, siz=siz)
	;oplot, rcs2
	rcs=(rcs+rcs2)/2.
	averages-=1
endwhile

;plot, rcs

deconvs = deconvolve_simple(data=rcs, iter=iter, og=og_source, matrix=matrix)

dc=size(deconvs)
deconvolutions = []
coors = []
for i=0, dc[2]-1 do begin
	d = deconvs[*,i]
	coors = [coors, d[-1]]
	deconvolutions = [[deconvolutions], [d[0:-2]]]
endfor


print, total(rcs[7*60:9*60])  
print, total(phi_psf[7*60:9*60])  
print, total(rcs[7*60:9*60])/total(phi_psf[7*60:9*60]) 

print, total(rcs)  
print, total(phi_psf)  
print, total(rcs)/total(phi_psf) 

;===============================================================================
;===============================================================================
;Now time to plot things
;===============================================================================
;===============================================================================

if keyword_set(FULLPLOT) then begin


	popen, 'model_'+source+'_'+method+'_'+erez+'_'+edist+'_'+'thresh'+strtrim(threshold, 2)+'.ps', $
			xsi=6, ysi=10
	!Y.margin=4.
	!X.margin=4.
	ch=1.1
	th=8
	lnth=4
	fth=4
	!p.multi=[0,1,5]

	th=5
	
	;===============================================================================
	;===============================================================================

	loadct, 3

	plot, og_source[indices], xtitle='X (Source Coordinates)', ytitle = 'Counts', thick=th, $
			title='Original Unconvolved Source', charsize=1.5, yrange=[0, max(og_source[indices])*1.2], $
			ytickinterval=50, xtickinterval=20, xthick=th, ythick=th, font=-1, charthick=fth
		
	oplot, og_source[indices], thick=th, color=150
		
	;===============================================================================

	toplot=	phi_psf[indices]
	
	plot, toplot, xtitle='X ('+strtrim(stripbins,2)+' Bins/Strip)', ytitle = 'Counts', thick=0, $
			title='Convolved Source (PSF)',linestyle=6, psym=0, charsize=1.5, $
			ytickinterval=50, xtickinterval=20, xthick=th, ythick=th, font=-1, charthick=fth, $
			yrange=[0, max(og_source[indices])*1.2]
	for n=0, n_elements(toplot)-1 do begin
		oplot, [n-0.5,n-0.5], [0, toplot[n]], thick=th
		oplot, [n-0.5,n-0.5]+1, [0, toplot[n]], thick=th
		oplot, [n-0.5, n+1-0.5], [toplot[n], toplot[n]], thick=th
	endfor
	n=0
	while n LT n_elements(phi_psf) do begin
	  oplot, [n,n,n], [0,1,4*max(phi_psf)], thick=th/4., linestyle=2
	  n+=stripbins
	endwhile

	loadct, 3

	oplot, og_source[indices], color=150, thick=th

	loadct, 8

	;===============================================================================

	toy = findgen(20)+0.5	

	toplot = sl_hist[indices_strip]
		;toy[indices_strip],
	plot, toplot, xtitle='X (Strips)', ytitle = 'Counts', thick=th, $
			title='Strip-binned Source (Post Charge Sharing)',linestyle=6, psym=4, charsize=1.5, $
			ytickinterval=max(toplot)/2, xtickinterval=1, symsize=0.1, xstyle=1, xthick=th, ythick=th, $
			font=-1, charthick=fth, yrange=[0, max(toplot)*1.5]
	for n=0, n_elements(toplot)-1 do begin
		oplot, [n,n], [0, toplot[n]], thick=th
		oplot, [n,n]+1, [0, toplot[n]], thick=th
		oplot, [n, n+1], [toplot[n], toplot[n]], thick=th
	endfor
	n=0
	while n LT n_elements(sl_hist) do begin
	  oplot, [n,n,n], [0,1,4*max(sl_hist)], thick=th/4., linestyle=2
	  n+=1
	endwhile

	;===============================================================================

	toplot=rcs[indices]

	plot, toplot, xtitle='X ('+strtrim(stripbins,2)+' Bins/Strip)', ytitle = 'Counts', thick=0, $
			title='Reconstructed Source',linestyle=6, psym=0, yrange=[0.,max(og_source[indices])*1.1], $
			charsize=1.5, xthick=th, ythick=th, font=-1, charthick=fth, $
			ytickinterval=50, xtickinterval=20
	for n=0, n_elements(toplot)-1 do begin
		oplot, [n-0.5,n-0.5], [0, toplot[n]], thick=th
		oplot, [n-0.5,n-0.5]+1, [0, toplot[n]], thick=th
		oplot, [n-0.5, n+1-0.5], [toplot[n], toplot[n]], thick=th
	endfor			
	n=0
	while n LT n_elements(phi_psf) do begin
	  oplot, [n,n,n], [0,1,4*max(phi_psf)], thick=th/4., linestyle=2
	  n+=stripbins
	endwhile

	loadct, 3

	oplot, og_source[indices], color=150, thick=th
	al_legend, ['Original Source', 'Reconstructed (pre-deconv)'], textcol=[150,0], /right, /top, /clear

	;===============================================================================

	loadct, 8

	mx_loc = where(coors EQ max(coors))

	samples=[0, 1, 2, 3]

	toplot= deconvolutions[*,0]

	plot, toplot[indices], ytitle = 'Counts', thick=th/2, $
			title='Deconvolved Reconstructed Source', yrange=[0.,max(og_source[indices])*1.1], $
			charsize=1.5, xtitle='X (Source Coordinates)', color=0, $
			ytickinterval=50, xtickinterval=20, xthick=th, ythick=th, font=-1, charthick=fth
	n=0

	oplot, rcs[indices], color=30, thick=4
	labels=['iter=0']
	colors=[30]

	loadct, 3

	oplot, og_source[indices], color=150, thick=th

	loadct, 8

	while n LT n_elements(toplot[indices]) do begin
	  oplot, [n,n,n], [0,1,4*max(phi_psf)], thick=th/4., linestyle=2
	  n+=stripbins
	endwhile
	color=60
	for i=0, n_elements(samples)-1 do begin
		ind = samples[i]
		colors=[colors, color]
		toplot= deconvolutions[*,ind]
		oplot, toplot[indices], color=color, thick=4
		color+=30
		;print, 'One Plot'
		labels=[labels, 'iter='+strtrim(string(round(ind+1)), 1)]
	endfor

	colors=[colors,0]
	labels=[labels, 'iter='+strtrim(string(round(mx_loc+1)), 1)+' ']

	mx_loc = where(coors EQ max(coors))
	toplot=deconvolutions[*,mx_loc[0]]
	oplot, toplot[indices], thick=8, color=0


	al_legend, labels, textcol=colors, box=1, /right, /top, /clear
	
	result = toplot[indices]
	ogsource = og_source[indices]
	save, result, ogsource, filename='results_for_evaluation.sav'

	;===============================================================================
	;===============================================================================

	!p.multi=0
	!Y.margin=[4.,2.]
	pclose
	spawn, 'open model_'+source+'_'+method+'_'+erez+'_'+edist+'_'+'thresh'+strtrim(threshold, 2)+'.ps'

endif



if keyword_set(COMPACTPLOT) then begin


	popen, 'compact_model_'+source+'_'+method+'_'+erez+'_'+edist+'_'+'thresh'+strtrim(threshold, 2)+'.ps', $
			xsi=6, ysi=10
			
	!Y.margin=4.
	!X.margin=4.
	ch=1.1
	th=8
	lnth=4
	fth=4
	!p.multi=[0,1,5]

	th=5

	;===============================================================================
	;===============================================================================

	loadct, 3

	plot, og_source[indices], xtitle='X (Source Coordinates)', ytitle = 'Counts', thick=th, $
			title='Original Unconvolved Source', charsize=1.5, yrange=[0, max(og_source[indices])*1.2], $
			ytickinterval=50, xtickinterval=20, xthick=th, ythick=th, font=-1, charthick=fth
		
	n=0	
	while n LT n_elements(og_source[indices]) do begin
	  oplot, [n,n,n], [0,1,4*max(og_source)], thick=th/4., linestyle=2
	  n+=stripbins
	endwhile
		
	oplot, og_source[indices], thick=th, color=150


		
	;===============================================================================
	;===============================================================================

	toy = findgen(20)+0.5	

	toplot = sl_hist[indices_strip]
		;toy[indices_strip],
	plot, toplot, xtitle='X (Strips)', ytitle = 'Counts', thick=th, $
			title='Strip-binned Source (Post Charge Sharing)',linestyle=6, psym=4, charsize=1.5, $
			ytickinterval=max(toplot)/2, xtickinterval=1, symsize=0.1, xstyle=1, xthick=th, ythick=th, $
			font=-1, charthick=fth, yrange=[0, max(toplot)*1.1];, xrange=[0,8]
	for n=0, n_elements(toplot)-1 do begin
		oplot, [n,n], [0, toplot[n]], thick=th
		oplot, [n,n]+1, [0, toplot[n]], thick=th
		oplot, [n, n+1], [toplot[n], toplot[n]], thick=th
	endfor
	n=0
	while n LT n_elements(sl_hist) do begin
	  oplot, [n,n,n], [0,1,4*max(sl_hist)], thick=th/4., linestyle=2
	  n+=1
	endwhile

	;===============================================================================
	;===============================================================================

	loadct, 8

	mx_loc = where(coors EQ max(coors))

	samples=[0, 1, 2, 3]

	toplot= deconvolutions[*,0]

	plot, toplot[indices], ytitle = 'Counts', thick=th/2, $
			title='Deconvolved Reconstructed Source', yrange=[0.,max(og_source[indices])*1.1], $
			charsize=1.5, xtitle='X (Source Coordinates)', color=0, $
			ytickinterval=50, xtickinterval=20, xthick=th, ythick=th, font=-1, charthick=fth
	n=0

	oplot, rcs[indices], color=30, thick=4
	labels=['iter=0']
	colors=[30]

	loadct, 3

	oplot, og_source[indices], color=150, thick=th

	loadct, 8

	while n LT n_elements(toplot[indices]) do begin
	  oplot, [n,n,n], [0,1,4*max(phi_psf)], thick=th/4., linestyle=2
	  n+=stripbins
	endwhile
	color=60
	for i=0, n_elements(samples)-1 do begin
		ind = samples[i]
		colors=[colors, color]
		toplot= deconvolutions[*,ind]
		oplot, toplot[indices], color=color, thick=4
		color+=30
		;print, 'One Plot'
		labels=[labels, 'iter='+strtrim(string(round(ind+1)), 1)]
	endfor

	colors=[colors,0]
	labels=[labels, 'iter='+strtrim(string(round(mx_loc+1)), 1)+' ']

	mx_loc = where(coors EQ max(coors))
	toplot=deconvolutions[*,mx_loc[0]]
	oplot, toplot[indices], thick=8, color=0


	al_legend, labels, textcol=colors, box=1, /right, /top, /clear
	
	
	

	;===============================================================================
	;===============================================================================
	;===============================================================================
	;===============================================================================


	!p.multi=0
	!Y.margin=[4.,2.]
	pclose
	spawn, 'open compact_model_'+source+'_'+method+'_'+erez+'_'+edist+'_'+'thresh'+strtrim(threshold, 2)+'.ps'
			
	
endif


if keyword_set(PLOTSPECTRA) then begin

	restore, 'saved_restored_spectrum.sav'
	restore, 'input_spectrum.sav'

	popen, 'spectrum_'+source+'_'+edist+'_'+'thresh_'+strtrim(threshold,2)+'.ps', xsi=7, ysi=10
	!Y.margin=4.
	!X.margin=4.
	ch=1.1
	th=6
	lnth=4
	!p.multi=[0,1,4]

	th=7
	fth=4
	;===============================================================================
	;===============================================================================

	binn=0.2
	plotextent = 25/binn

	is = histogram(input_spectrum, binsize=binn, nbins=plotextent, min=0, locations=is_x)
	rs = histogram(raw_spectrum, binsize=binn, nbins=plotextent, min=0, locations=rs_x)
	ss = histogram(spectrum, binsize=binn, nbins=plotextent, min=0, locations=ss_x)
	cs = histogram(rawcsa_spectrum, binsize=binn, nbins=plotextent, min=0, locations=cs_x)
	cps = histogram(rec_comp_spectrum, binsize=binn, nbins=plotextent, min=0, locations=cps_x)

	plot, is_x, is, title='Input Spectrum', linestyle=6, psym=4, $
			ytitle='Counts/bin', xthick=th, ythick=th, font=-1, charthick=fth, charsize=2, $
			yrange=[0,1.2*max(is)], ystyle=1, symsize=0.1, xtitle='Energy (keV)', xtickinterval=2
	for n=0, n_elements(is)-1 do begin
		oplot, [is_x[n],is_x[n]], [0, is[n]], thick=th
		oplot, [is_x[n],is_x[n]]+binn, [0, is[n]], thick=th
		oplot, [is_x[n], is_x[n]+binn], [is[n], is[n]], thick=th
	endfor
	oplot, [threshold, threshold], [0, 10000], color=50, thick=th
	oplot, [sec_threshold, sec_threshold], [0, 10000], color=50, thick=th, linestyle=2

	if edist EQ 'ALS' then begin
		oplot, [7,7], [0,10000], color=100, thick=th
		oplot, [21,21], [0,10000], color=100, thick=th
		oplot, [28,28], [0,10000], color=100, thick=th
	endif

	if edist EQ 'ALS_high' then begin
		oplot, [21,21], [0,10000], color=100, thick=th
	endif


	plot, rs_x, rs, title='Charge Shared Spectrum', linestyle=6, psym=4, $
			ytitle='Counts/bin', xthick=th, ythick=th, font=-1, charthick=fth, charsize=2, $
			yrange=[0,1.2*max(rs)], ystyle=1, symsize=0.1, xtitle='Energy (keV)', xtickinterval=2
	for n=0, n_elements(rs)-1 do begin
		oplot, [rs_x[n],rs_x[n]], [0, rs[n]], thick=th
		oplot, [rs_x[n],rs_x[n]]+binn, [0, rs[n]], thick=th
		oplot, [rs_x[n], rs_x[n]+binn], [rs[n], rs[n]], thick=th
	endfor
	oplot, [threshold, threshold], [0, 10000], color=50, thick=th
	oplot, [sec_threshold, sec_threshold], [0, 10000], color=50, thick=th, linestyle=2

	if edist EQ 'ALS' then begin
		oplot, [7,7], [0,10000], color=100, thick=th
		oplot, [21,21], [0,10000], color=100, thick=th
		oplot, [28,28], [0,10000], color=100, thick=th
	endif
	
	if edist EQ 'ALS_high' then begin
		oplot, [21,21], [0,10000], color=100, thick=th
	endif

	plot, cs_x, cs, title='Reconstructed Spectrum', linestyle=6, psym=4, $
			ytitle='Counts/bin', xthick=th, ythick=th, font=-1, charthick=fth, charsize=1.75, $
			yrange=[0,1.2*max(cs)], ystyle=1, symsize=0.1, xtitle='Energy (keV)', xtickinterval=2
	for n=0, n_elements(cs)-1 do begin
		oplot, [cs_x[n],cs_x[n]], [0, cs[n]], thick=th
		oplot, [cs_x[n],cs_x[n]]+binn, [0, cs[n]], thick=th
		oplot, [cs_x[n], cs_x[n]+binn], [cs[n], cs[n]], thick=th
	endfor
	oplot, [threshold, threshold], [0, 10000], color=50, thick=th
	oplot, [sec_threshold, sec_threshold], [0, 10000], color=50, thick=th, linestyle=2
	
	if edist EQ 'ALS' then begin
		oplot, [7,7], [0,10000], color=100, thick=th
		oplot, [21,21], [0,10000], color=100, thick=th
		oplot, [28,28], [0,10000], color=100, thick=th
	endif

	if edist EQ 'ALS_high' then begin
		oplot, [21,21], [0,10000], color=100, thick=th
	endif

	;===============================================================================
	;===============================================================================

	!p.multi=0
	!Y.margin=[4.,2.]
	pclose
	spawn, 'open spectrum_'+source+'_'+edist+'_'+'thresh_'+strtrim(threshold,2)+'.ps'


	;===============================================================================
	;===============================================================================



endif


stop

end




















