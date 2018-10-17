World Magnetic Model WMM2015v2
========================================================
Date September 24, 2018

WMM.COF					WMM2015v2 Coefficients file
						(Replace old WMM.COF (WMM2015 
						 or WMM2010) file with this)

1. Installation Instructions
==========================

WMM2015 GUI
-----------

Go to installed directory, find WMM.COF and remove it. 
Replace it with the new WMM.COF. 


WMM_Linux and WMM_Windows C software
------------------------------------

Replace the WMM.COF file in the "bin" directory with the
new provided WMM.COF file 


Your own software
-----------------
Depending on your installation, find the coefficient 
file, WMM.COF (unless renamed). Replace it with the new
WMM.COF file (renaming it appropriately if necessary).

If the coefficients are embedded in your software, you
may need to embed the new coefficients.

2. Installation Verification
============================

To confirm you are using the correct WMM.COF file open 
it and verify that the header is:

    2015.0            WMM-2015v2      09/18/2018

To assist in confirming that the installation of the new 
coefficient file is correct we provide a set of test 
values in this package.  Here are a few as an example:

Date     HAE   Lat   Long   Decl   Incl       H         X         Y        Z         F        Ddot   Idot   Hdot   Xdot   Ydot   Zdot   Fdot
2015.0   100     0    120    0.54 -15.98   37539.7   37538.1     351.1  -10751.1   39048.9   -0.07   0.09   18.1   18.5  -46.1   60.2    0.8
2015.0   100   -80    240   69.22 -72.57   15818.3    5612.2   14789.3  -50385.8   52810.5   -0.08   0.04   10.3   25.1    1.5   82.5  -75.6
2017.5     0    80      0   -2.59  83.08    6612.0    6605.2    -298.7   54506.3   54905.9    0.53   0.02  -15.3  -12.6   61.3   39.0   36.9
2017.5     0     0    120    0.37 -15.63   39570.2   39569.4     252.3  -11067.9   41088.9   -0.08   0.09   19.0   19.3  -50.2   64.3    1.0

Where HAE is height above WGS-84 ellipsoid.

To confirm you are using the correct test value file 
please check the title, which should read 
"Test Values for WMM2015v2".



Model Software Support
======================

*  National Centers for Environmental Information (NCEI)
*  E/NE42 325 Broadway
*  Boulder, CO 80305 USA
*  Attn: Manoj Nair or Arnaud Chulliat
*  Phone:  (303) 497-4642 or -6522
*  Email:  geomag.models@noaa.gov
For more details about the World Magnetic Model visit 
http://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml


