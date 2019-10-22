;+
; NAME:
;   fit1d_acs_profile_nuker_template.pro
;
; AUTHOR:
;   Arna Karick, Astrophysics Group, University of Oxford, UK.
;
; CONTACT: See www.arnakarick.com (Twitter @drarnakarick)
;
; ORIGINAL VERSION: Karick, A.M. and Carter, D. (HST Coma Cluster
;   Treasury Survey Team), June 2011. 
;
; THIS VERSION: This version created June 2012 and used in "ATLAS3D Project
;   – XXIII. Angular momentum and nuclear surface brightness profiles."
;   (Krajnovic, D., Karick, A.M., Davies, R.L, Thorsten, N., Sarzi, M,
;   Emsellem, E., Cappellari, M., et al. 2013,MNRAS, 433, 2812, doi:
;   10.1093/mnras/stt905)
;  
; 
; Uses MPFITFUN.PRO: Levenberg-Marquardt least-squares fit to IDL functions
;
;
; MPFITFUN AUTHOR
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; USE:
;   The user must generate isophote tables using the
;   IRAF.stsdas.analysis.isophote.ellipse package.
;
; RUN: standard IDL command:
; idl> .run fit1d_acs_profile_nuker.pro
; idl> fit1d_acs_profile_nuker
;

pro fit1d_acs_profile_nuker,ps=ps
;,ps=ps

close, /all & erase

; --------------------------------- read in Atlas3D profiles and print ---------------------
;                                   calibrated parameters to a new file 
fileA='NGCNUM_FILT.dat'
  readcol,fileA, $
  sma,intens,intens_err,ellip,ellip_err,pa,pa_err,a3,a3_err,b3,b3_err,a4,a4_err,b4,b4_err, $
      format='(F,F,F,X,X,F,F,F,F,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,X,F,F,F,F,F,F,F,F)',$
      skipline=4

fileB='NGCNUM_FILT_conv.dat'
  readcol,fileB, $
  csma,cintens,cintens_err, $
      format='(F,F,F,F)', $
      skipline=4

; calculate geometric mean

sma_asec = sma * 0.05
area_pix = 0.05 ^ 2.0
;rgeo = sma_asec * sqrt(1-ellip)  ; rgeo = sqrt (rgeo)
rgeo = sma_asec   ; use the sma - but call it rgeo in the script

;ZPTMAG_F475W = 26.04541
ZPTMAG_F475W = 46.87 - 21.10 - 5.0*alog10(4745) + 18.6921
ZPTMAG_F555W = 46.79 - 21.10 - 5.0*alog10(5360) + 18.6921
ZPTMAG_F606W = 44.77 - 21.10 - 5.0*alog10(5919) + 18.6921
ZPTMAG_F814W = 47.89 - 21.10 - 5.0*alog10(8059) + 18.6921

;intens_sb  = -2.5 * alog10(intens[i]) + 26.04541 + VZPT - 2.5 *
;                                                          alog10(1/area_pix) ;vzpt color term
intens_sb  = -2.5*alog10(intens)  - 2.5*alog10(1/area_pix) + ZPTMAG_FILT 
intens_sb_err =  -2.5 * 0.434 * (intens_err/intens)

cintens_sb = -2.5*alog10(cintens) - 2.5*alog10(1/area_pix) + ZPTMAG_FILT
cintens_sb_err =  -2.5 * 0.434 * (cintens_err/cintens)

decon_sb = 2*intens_sb - cintens_sb
decon_sb_err = sqrt(intens_sb_err^2 + cintens_sb_err^2)


;par=[muRb,a,g,b,Rb]
;par=[Ib,a,g,b,Rb]

pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},5)

; ---------------------- limits for mu(Rb)
pi(0).limited(0)=1
pi(0).limited(1)=1
pi(0).limits(0) = 14.D
pi(0).limits(1) = 24.D

; ---------------------- limits for alpha
pi(1).limited(0)=1
pi(1).limited(1)=1
pi(1).limits(0) = 0.D
pi(1).limits(1) = 5.D
;pi(1).fixed=1
;initial(1)=5.0
;
; ---------------------- limits for gamma - inner power-law slope
pi(2).limited(0)=1
pi(2).limited(1)=1
pi(2).limits(0) = -3.D
pi(2).limits(1) = 5.D
;pi(2).fixed=1
;initial(2)=0.081

; ---------------------- limits for beta - outer power-law slope
pi(3).limited(0)=1
pi(3).limited(1)=1
pi(3).limits(0) = 0.D
pi(3).limits(1) = 10.D
;pi(3).fixed=1

; ---------------------- limits for Rb - break radius
pi(4).limited(0)=1  
pi(4).limited(1)=1  
pi(4).limits(0) = 0.1D
pi(4).limits(1) = 100.D
;pi(4).fixed=1
;start(4)=1.46

start=[17.0,1.0,1.0,5.0,2.0]

rmin = RMIN
rmax = RMAX

psf = where(rgeo[*] gt rmin and rgeo[*] lt rmax)

result=mpfitfun('nuker',rgeo[psf],decon_sb[psf],decon_sb_err[psf],start,STATUS=status,ERRMSG=errmsg,weights=1D,parinfo=pi,perror=perror)  ;  eual weights

if status LE 0 then message, errmsg

;fit_Id = result(0)
fit_mub = result(0)
fit_a  = result(1)
fit_g  = result(2)
fit_b  = result(3)
fit_Rb = result(4)

print,'mu(Ib) (at core.break radius)=',result(0)
print,'alpha)=',result(1)
print,'gamma (inner power-law slope)=',result(2)  ; g<0.3 "core galaxies", g>0.5 "power-law" galaxies (Lauer et al.1995)
print,'beta (outer power-law slope)=',result(3)
print,'Rb (break radius)=',result(4)

; Convert Rb break radius into parsecs
dist = DMPC
fit_Rb_pc = (fit_Rb/3600.0)*(!pi/180.0)*dist*1000000.0
res_pc    = (0.049/3600.0)*(!pi/180.0)*dist*1000000.0
print,'Rb (break radius in parsecs)=',fit_Rb_pc
print,'image/pixel resolution in parsecs=',res_pc
rmin_pc = (rmin/3600.0)*(!pi/180.0)*dist*1000000.0
rmax_pc = (rmax/3600.0)*(!pi/180.0)*dist*1000000.0

fit = fit_mub - 2.5*(fit_b-fit_g)/fit_a*alog10(2) - 2.5*fit_g*alog10(fit_Rb/rgeo) - 2.5*(fit_g-fit_b)/fit_a*alog10(1+(rgeo/fit_Rb)^fit_a)

diff = fit - decon_sb

; calculate the RMS off the difference between data and fit
diff_rlim = diff(where(rgeo[*] gt rmin and rgeo[*] lt rmax))
diff_rms  = sqrt(total(diff_rlim^2)/n_elements(diff_rlim))

rdash = RDASH
gamma_dash = (fit_g + fit_b*((rdash/fit_Rb)^fit_a))/(1+(rdash/fit_Rb)^fit_a)

print,'gamma_dash (local gradient at',rdash,' ")',gamma_dash

;if (keyword_set(ps)) then $
; apj_open,'fig_'+ngc[i]+'_1dS',/encap,/color,xsize=4.0,/inches,ysi=4.0

loadct,2
y1=min(decon_sb)
y2=max(decon_sb)
ploterror,rgeo,decon_sb,decon_sb_err,psym=4, yrange=[y2,y1], /xlog
oplot,rgeo,fit, color=30, thick=2
oplot,rgeo[psf],fit[psf], color=60, thick=2


openw, unit2,repstr(fileA,'.dat','_decon_1dN.dat'), width=200,bufsize=0,/get_lun

printf, unit2, '# fit_Rb(pc) fit_Rb    fit_mub     fit_a      fit_b     fit_g  gamma_dash      rmin      rmin_pc  rmax   rmax_pc   rms      res_pc',format='(a130)'  
printf, unit2, fit_Rb_pc, fit_Rb, fit_mub, fit_a,  fit_b, fit_g, gamma_dash, rmin, rmin_pc, rmax, rmax_pc, diff_rms, res_pc, format='(f8.2,3x,f6.2,3x,f8.2,3x,f8.2,3x,f8.2,3x,f8.2,3x,f8.2,3x,f8.2,3x,f8.2,3x,f5.2,3x,f8.2,3x,f6.2,3x,f6.2,3x,f6.2)'
printf, unit2, '#    sma(asec)      rgeo          decon_sb        decon_sb_err  intens_sb      intens_sb_err   ellip          ellip_err     pa              pa_err         fit                       diff',format='(a186)'

nfitpoints = n_elements(fit)

for j=0,nfitpoints-1 do begin

;printf, unit2, sma[psf][j],rgeo[psf][j],sb[psf][j],fit[psf][j],diff[psf][j],ellip[psf][j],pa[psf][j],format='(F,F,F,F,F,F,F)'
printf, unit2,sma_asec[j],rgeo[j],decon_sb[j],decon_sb_err[j],intens_sb[j],intens_sb_err[j],ellip[j],ellip_err[j],pa[j],pa_err[j],fit[j],diff[j],format='(F,F,F,F,F,F,F,F,F,F,F,F)'

endfor

free_lun,unit2
close, /all

end
