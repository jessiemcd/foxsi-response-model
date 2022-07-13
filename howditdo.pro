pro howditdo


;Evaluating the success of the reconstruction process (post deconvolution - output from response_modeling)
;Somewhat hardcoded to examples from thesis

restore, 'results_for_evaluation_edges_2thresh.sav'
source1=ogsource - make_array(n_elements(ogsource), val=15.) 
result1=result - make_array(n_elements(ogsource), val=15.) 

restore, 'results_for_evaluation_centers_2thresh.sav'
source2=ogsource - make_array(n_elements(ogsource), val=15.) 
result2=result - make_array(n_elements(ogsource), val=15.) 

stripbins=60
strips=20
siz = strips*stripbins
xis = indgen(siz)

edge=0
;Generating observed source: two peaks on continuum
loc1 = 7.*stripbins
loc2 = 10.*stripbins
if edge EQ 0 then BEGIN
	loc1+=0.5*stripbins
	loc2+=0.5*stripbins
endif
sigma1 = 1./10*stripbins
sigma2 = 2./10*stripbins
mag1 = 50.
mag2 = 70.

vall=0.
psi_xi = mag1*EXP(-((xis-loc1)^2/(2.*sigma1^2))) + $
	  mag2*EXP(-((xis-loc2)^2/(2.*sigma2^2)))
psi_xi = psi_xi + make_array(n_elements(xis), val=vall)  


print, 'Center 1: ', loc1, ' Center 2:', loc2
print, 'Width 1: ', sigma1, ' Width 2:', sigma2
print, 'Width 1 FWHM: ', sigma1*2.355, ' Width 2 FWHM:', sigma2*2.355

plot, psi_xi, thick=th

peak=where(psi_xi EQ max(psi_xi))
peak=peak[0]
iv=20.

xx=xis

;make data interval for fitting first peak + do gaussian fit
fitbins = xx[where(xx GE peak-iv)]
fitbins2 = fitbins[where(fitbins LE peak+iv)]
fithist = psi_xi[where(xx GE peak-iv)]
fithist2 = fithist[where(fitbins LE peak+iv)]
fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

;plot fit data interval, fit result, and fit center 
oplot, fitbins2, fithist, thick=th, col=180,  psym=10
oplot, fitbins2, fit, thick=th, col=180, linestyle=2
oplot, [coeff[1],coeff[1]], [0.1, max(psi_xi)*2], linestyle=2, col=100, thick=th

print, 'Source Peak Amplitude: ', strtrim(coeff[0],2), ' Center: ', strtrim(coeff[1],2), ' Sigma: ', strtrim(coeff[2],2)





popen, 'howditdo.ps', $
		xsi=8, ysi=10
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

loadct, 2

xx=indgen(n_elements(result1))

r1maxes = lclxtrem(result1, 5, /maxima)
print, 'Result Edges Maxes: ', r1maxes
s1maxes = lclxtrem(source1, 5, /maxima)
print, 'Source Edges Maxes: ', s1maxes

r2maxes = lclxtrem(result2, 5, /maxima)
print, 'Result Centers Maxes: ',r2maxes
s2maxes = lclxtrem(source2, 5, /maxima)
print, 'Source Centers Maxes: ',s2maxes

plot, result1, thick=th, yrange=[0, max(result1)*1.5], title='Edges Example'
oplot, source1, color=30, thick=th

;for i=0, n_elements(r1maxes)-1 do begin
;	oplot, [r1maxes[i], r1maxes[i]], [0,max(result1)*1.5]
;endfor
;for i=0, n_elements(s1maxes)-1 do begin
;	oplot, [s1maxes[i], s1maxes[i]], [0,max(result1)*1.5], color=30
;endfor

sourcepeaks=[]
resultpeaks=[]
sourcewidths=[]
resultwidths=[]
resultamps = []
sourceamps = []

for i=0, 1 do begin

	peak=s1maxes[i]
	iv=20

	;make data interval for fitting first peak + do gaussian fit
	fitbins = xx[where(xx GE peak-iv)]
	fitbins2 = fitbins[where(fitbins LE peak+iv)]
	fithist = source1[where(xx GE peak-iv)]
	fithist2 = fithist[where(fitbins LE peak+iv)]
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

	;plot fit data interval, fit result, and fit center 
	oplot, fitbins2, fithist, thick=th, col=180,  psym=10
	oplot, fitbins2, fit, thick=th, col=180, linestyle=2
	oplot, [coeff[1],coeff[1]], [0.1, max(source1)*2], linestyle=2, col=100, thick=th
	
	print, 'Source Peak '+strtrim(i, 2)+':', ' Amplitude: ', strtrim(coeff[0],2), ' Center: ', strtrim(coeff[1],2), ' Sigma: ', strtrim(coeff[2],2)
	
	sourcepeaks=[sourcepeaks, coeff[1]]
	sourcewidths=[sourcewidths, coeff[2]]
	sourceamps = [sourceamps, coeff[0]]
	
endfor

print, ''

for i=0, 3 do begin

	peak=r1maxes[i]
	iv=15

	;make data interval for fitting first peak + do gaussian fit
	fitbins = xx[where(xx GE peak-iv)]
	fitbins2 = fitbins[where(fitbins LE peak+iv)]
	fithist = result1[where(xx GE peak-iv)]
	fithist2 = fithist[where(fitbins LE peak+iv)]
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

	;plot fit data interval, fit result, and fit center 
	oplot, fitbins2, fithist, thick=th, col=150,  psym=10
	oplot, fitbins2, fit, thick=th, col=150, linestyle=2
	oplot, [coeff[1],coeff[1]], [0.1, max(source1)*2], linestyle=2, col=150, thick=th
	
	print, 'Result Peak '+strtrim(i, 2)+':', ' Amplitude: ', strtrim(coeff[0],2), ' Center: ', strtrim(coeff[1],2), ' Sigma: ', strtrim(coeff[2],2)
	
	resultpeaks=[resultpeaks, coeff[1]]
	resultwidths=[resultwidths, coeff[2]]
	resultamps = [resultamps, coeff[0]]
	
endfor

print, 'Peaks: '
print, sourcepeaks
print, resultpeaks
print, sourcepeaks-resultpeaks[0:1]

print, 'Widths: '
print, sourcewidths*2.355
print, resultwidths*2.355
print, (resultwidths[0:1]/sourcewidths)
;print, sourcewidths
;print, resultwidths
;print, (sourcewidths-resultwidths)


print, 'Amplitudes: '
print, sourceamps*2.355
print, resultamps*2.355
print, (resultamps[0:1]/sourceamps)
print, 'Highest Artifact Amplitude: ', resultamps[2]/resultamps[1]


print, ''
print, ''

sourcepeaks=[]
resultpeaks=[]
sourcewidths=[]
resultwidths=[]
resultamps = []
sourceamps = []

plot, result2, thick=th, yrange=[0, max(result2)*1.5], title='Centers Example'
oplot, source2, color=30, thick=th

;for i=0, n_elements(r2maxes)-1 do begin
;	oplot, [r2maxes[i], r2maxes[i]], [0,max(result2)*1.5]
;endfor
;for i=0, n_elements(s2maxes)-1 do begin
;	oplot, [s2maxes[i], s2maxes[i]], [0,max(result2)*1.5], color=30
;endfor

for i=0, 1 do begin

	peak=s2maxes[i]
	iv=20

	;make data interval for fitting first peak + do gaussian fit
	fitbins = xx[where(xx GE peak-iv)]
	fitbins2 = fitbins[where(fitbins LE peak+iv)]
	fithist = source2[where(xx GE peak-iv)]
	fithist2 = fithist[where(fitbins LE peak+iv)]
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

	;plot fit data interval, fit result, and fit center 
	oplot, fitbins2, fithist, thick=th, col=180,  psym=10
	oplot, fitbins2, fit, thick=th, col=180, linestyle=2
	oplot, [coeff[1],coeff[1]], [0.1, max(source2)*2], linestyle=2, col=100, thick=th
	
	print, 'Source Peak '+strtrim(i, 2)+':', ' Amplitude: ', strtrim(coeff[0],2), ' Center: ', strtrim(coeff[1],2), ' Sigma: ', strtrim(coeff[2],2)
	
	sourcepeaks=[sourcepeaks, coeff[1]]
	sourcewidths=[sourcewidths, coeff[2]]
	sourceamps = [sourceamps, coeff[0]]
	
endfor

print, ''

for i=0, 3 do begin

	peak=r2maxes[i]
	iv=15

	;make data interval for fitting first peak + do gaussian fit
	fitbins = xx[where(xx GE peak-iv)]
	fitbins2 = fitbins[where(fitbins LE peak+iv)]
	fithist = result2[where(xx GE peak-iv)]
	fithist2 = fithist[where(fitbins LE peak+iv)]
	fit=GAUSSFIT(fitbins2,fithist2,coeff,nterms=3)

	;plot fit data interval, fit result, and fit center 
	oplot, fitbins2, fithist, thick=th, col=150,  psym=10
	oplot, fitbins2, fit, thick=th, col=150, linestyle=2
	oplot, [coeff[1],coeff[1]], [0.1, max(source2)*2], linestyle=2, col=150, thick=th
	
	print, 'Result Peak '+strtrim(i, 2)+':', ' Amplitude: ', strtrim(coeff[0],2), ' Center: ', strtrim(coeff[1],2), ' Sigma: ', strtrim(coeff[2],2)
	
	resultpeaks=[resultpeaks, coeff[1]]
	resultwidths=[resultwidths, coeff[2]]
	resultamps = [resultamps, coeff[0]]
	
endfor

print, 'Peaks: '
print, sourcepeaks
print, resultpeaks
print, sourcepeaks-resultpeaks[0:1]

print, 'Widths: '
print, sourcewidths*2.355
print, resultwidths*2.355
print, (resultwidths[0:1]/sourcewidths)

print, 'Amplitudes: '
print, sourceamps*2.355
print, resultamps*2.355
print, (resultamps[0:1]/sourceamps)

print, 'Highest Artifact Amplitude: ', resultamps[2]/resultamps[1]

;===============================================================================
;===============================================================================

!p.multi=0
!Y.margin=[4.,2.]
pclose
spawn, 'open howditdo.ps'











stop
end