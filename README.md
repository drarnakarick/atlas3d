# The ATLAS3D Project – XXII

This directory contains sample data, and code from my research at University of Oxford.

The published journal paper:https://arxiv.org/abs/1305.4973

### IDL routine for fitting Nuker Profiles

This directory contains an example IDL routine for fitting Nuker
profiles to Hubble Space Telescope (ACS) images of Atlas3D
galaxies. The origianl routines were developed for the HST/ACS Coma
Cluster Treasury Survey and modified by the original author (Arna
Karick, June 2012) to analyse the HST galaxy images presented in
"ATLAS3D Project – XXIII. Angular momentum and nuclear surface
brightness profiles." (Krajnovic, D., Karick, A.M., Davies, R.L,
Thorsten, N., Sarzi, M, Emsellem, E., Cappellari, M., et al. 2013,
MNRAS, 433, 2812, doi: 10.1093/mnras/stt905)

We would appreciate it if researchers using this script, or
modifications of the script, acknowledge the original authors.

VERSION:
Original version: Karick, A.M. and Carter, D. (HST Coma Cluster
Treasury Survey Team), June 2011.  

This modified version: Karick, A.M., June 2012.

OPTIONAL: 
- IRAF: STSCI/Gemini Ureka installation (recommended); http://ssb.stsci.edu/ureka/

REQUIREMENTS:
- IDL:(requires a licence) http://www.exelisvis.com/ProductsServices/IDL.aspx  
- IDLAstro: http://idlastro.gsfc.nasa.gov


**1. ISOPHOTE FITTING (IRAF)**

The following is a brief summary of the method used to generate galaxy
isophotes (refer to 'Section 3.1 Surface brightness profiles' in
Krajnovic, D., Karick, et al. (2013). The
IRAF.stsdas.analysis.isophtoe.ellipse package was used manually to
generate isphotes.

1. Galaxy centre determined, minimum isophote radii set. Position
   angle and ellipticity allowed to vary throughout the process (IRAF
   ellipse and controlpar).
2. First pass isophotes created (.tab files) and overlaid on images as
   a quality control check (IRAF ispimap).
3. First pass isophotes used to generate a model image of the galaxy
   (IRAF bmodel).
4. Residual images created by subtracting model image (IRAF
   imarith). This was done to identify large features, identify bright
   stars, star forming reasons and dust features that might affect the
   isophote fitting.
5. Object masks created using the residual image to mask out spurious features (IRAF nproto.objmasks).
6. Masks edited to remove any flagged pixels at the centre of the image (IRAF imreplace).
7. First pass isophotes deleted and new ones created using the object masks. Output: NGCNUM_FILT.tab
8. Dealing with the image PSF: Create NGCNUM_FILT_conv.fits images by
   convolving the original image with a PSF model image - created
   using TinyTIM (Krist & Hook 2011) -
   http://www.stsci.edu/hst/observatory/focus/TinyTim
9. Create isophotes for this images using the previous isphote table (NGCNUM_FILT.tab) and mask files. 
10. PSF de-convolution and filter calibration is dealt with in the IDL script.

**2. FITTING ROUTINE (IDL)**

This routine uses the IDL routine mpfitfun.pro:
http://www.physics.wisc.edu/~craigm/idl/down/mpfitfun.pro by Craig
B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770,
craigm@lheamail.gsfc.nasa.gov UPDATED VERSIONs can be found on my WEB
PAGE: http://cow.physics.wisc.edu/~craigm/idl/idl.html

Output files form the ISOPHOTE fitting (converted to .dat files using
IRAF.ttools) file. This is the input data file for the IDL fitting
script.

HST Zerpoint magnitudes were determined following Sirianni et
al. (2005) and HST/ACS flux calibration zerpoint keywords:
http://www.stsci.edu/hst/acs/analysis/zeropoints

To explore the core region of galaxies we fit the ‘Nuker law’ (Lauer
et al. 1995), to our surface brightness profiles. For a detailed
discussion about the fitting equations used, refer to Sections 3.1 and
3.2 of Krajnovic, D., Karick, et al. (2013)

