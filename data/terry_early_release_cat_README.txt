This is an early-version high-level science product from GO-17776,
"A Precursor Survey of the Roman Galactic Bulge Time Domain Fields" 
for 8 individual fields in the Galactic bulge (HD16-ACS, HD16-WFC3UV,
HD70-ACS, HD70-WFC3UV, HD98-ACS, HD98-WFC3UV, HD138-ACS, and HD138-WFC3UV).
We provide ACS and WFC3 photometry, astrometry, as well as photometry 
for cross-identified sources in GO-17923 (WFC3-IR) and Gaia DR3. The 
photometry for every detected source is provided in at least 2 bands: 
F555W (V), F814W (I) using the VEGA mag photometric system. Where available, 
photometry for cross-identified stars is provided in near-IR filters F098M, 
F139M, F153M, F167N, and F130N. The astrometry is given in distortion-free pixel 
coordinates and equatorial coordinates (degrees, J2000). Drizzled 
images for this program can be accessed via MAST as part of Hubble 
Advanced Products (HAP). This early-release dataset is assigned the 
DOI: 10.17909/fq12-f295.

An overview of the program can be found in:
Terry et al. 2026 (submitted).
*******************************************************
Files:

terry_hst_bulge_survey_early_fields_cat.fits - FITS table of all 
detected sources in the eight early-release fields. A total of 767,024 
sources are in the table.

*******************************************************

The early release FITS table contains the following columns:

HST_ID = assigned ID of source

X = x position in the image (pixels)
X_err = error on x position in the image (pixels)

Y = y position in the image (pixels)
Y_err = error on y position in the image (pixels)                                     

F606W_Vegamag = calibrated magnitude in F606W (V)
F606W_mag_err = error on magnitude in F606W (V)

F814W_Vegamag = calibrated magnitude in F814W (I)
F814W_mag_err = error on magnitude in F814W (I)

N_detect = number of exposures detected in (max = 4)

HST_Field = Field name and detector (*_ACS or *_WFC3UV)

F098M_Vegamag = calibrated magnitude in F098M (WFC3-IR)
F098M_mag_err = error on magnitude in F098M

F139M_Vegamag = calibrated magnitude in F139M (WFC3-IR)
F139M_mag_err = error on magnitude in F139M

F153M_Vegamag = calibrated magnitude in F153M (WFC3-IR)
F153M_mag_err = error on magnitude in F153M

F167N_Vegamag = calibrated magnitude in F167N (WFC3-IR)
F167N_mag_err = error on magnitude in F167N

F130N_Vegamag = calibrated magnitude in F130N (WFC3-IR)
F130N_mag_err = error on magnitude in F130N

GaiaDR3_ID = ID of source in Gaia DR3

Gmag = Gaia G magnitude

BPmag = Gaia BP magnitude

RPmag = Gaia RP magnitude

RA = right ascension (J2000, degrees)
                                          
DEC = declination (J2000, degrees) 

*******************************************************

Note that the WFC3-IR photometry was reduced using DOLPHOT. 
The full reduction and calibration method for the WFC3-IR data 
will be presented in Nataf et al. (in prep).
