/*

   calin/math/constants.hpp -- Stephen Fegan -- 2015-11-05

   Namespace for constant values.

   Copyright 2015, Stephen Fegan <sfegan@llr.in2p3.fr>
   Laboratoire Leprince-Ringuet, CNRS/IN2P3, Ecole Polytechnique, Institut Polytechnique de Paris

   This file is part of "calin"

   "calin" is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 or
   later, as published by the Free Software Foundation.

   "calin" is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

*/

/*! \file Constants.hpp
  relativistic particle class header file

  \author   Stephen Fegan             \n
            UCLA                      \n
	    sfegan@astro.ucla.edu     \n

  \date     11/13/2004
  \version  1.0
  \note
*/

#pragma once

#include<iostream>
#include<iomanip>
#include<sstream>
#include<vector>
#include<set>
#include<string>
#include<limits>
#include<cmath>

namespace calin { namespace math { namespace constants {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Table of physical constants from NIST, included for reference as a comment
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*

  Fundamental Physical Constants --- Complete Listing

  From:  http://physics.nist.gov/constants

  Source: Peter J. Mohr and Barry N. Taylor, CODATA Recommended Values of the
          Fundamental Physical Constants: 2002,
          To be published.


  Quantity                                               Value                 Uncertainty          Unit
------------------------------------------------------------------------------------------------------------------------
speed of light in vacuum                               299 792 458           (exact)               m s^-1
magn. constant                                         12.566 370 614...e-7  (exact)               N A^-2
electric constant                                      8.854 187 817...e-12  (exact)               F m^-1
characteristic impedance of vacuum                     376.730 313 461...    (exact)               ohm
Newtonian constant of gravitation                      6.6742e-11            0.0010e-11            m^3 kg^-1 s^-2
Newtonian constant of gravitation over h-bar c         6.7087e-39            0.0010e-39            (GeV/c^2)^-2
Planck constant                                        6.626 0693e-34        0.000 0011e-34        J s
Planck constant in eV s                                4.135 667 43e-15      0.000 000 35e-15      eV s
Planck constant over 2 pi times c in MeV fm            197.326 968           0.000 017             MeV fm
Planck constant over 2 pi                              1.054 571 68e-34      0.000 000 18e-34      J s
Planck constant over 2 pi in eV s                      6.582 119 15e-16      0.000 000 56e-16      eV s
Planck mass                                            2.176 45e-8           0.000 16e-8           kg
Planck temperature                                     1.416 79e32           0.000 11e32           K
Planck length                                          1.616 24e-35          0.000 12e-35          m
Planck time                                            5.391 21e-44          0.000 40e-44          s
elementary charge                                      1.602 176 53e-19      0.000 000 14e-19      C
elementary charge over h                               2.417 989 40e14       0.000 000 21e14       A J^-1
magn. flux quantum                                     2.067 833 72e-15      0.000 000 18e-15      Wb
conductance quantum                                    7.748 091 733e-5      0.000 000 026e-5      S
inverse of conductance quantum                         12 906.403 725        0.000 043             ohm
Josephson constant                                     483 597.879e9         0.041e9               Hz V^-1
von Klitzing constant                                  25 812.807 449        0.000 086             ohm
Bohr magneton                                          927.400 949e-26       0.000 080e-26         J T^-1
Bohr magneton in eV/T                                  5.788 381 804e-5      0.000 000 039e-5      eV T^-1
Bohr magneton in Hz/T                                  13.996 2458e9         0.000 0012e9          Hz T^-1
Bohr magneton in inverse meters per tesla              46.686 4507           0.000 0040            m^-1 T^-1
Bohr magneton in K/T                                   0.671 7131            0.000 0012            K T^-1
nuclear magneton                                       5.050 783 43e-27      0.000 000 43e-27      J T^-1
nuclear magneton in eV/T                               3.152 451 259e-8      0.000 000 021e-8      eV T^-1
nuclear magneton in MHz/T                              7.622 593 71          0.000 000 65          MHz T^-1
nuclear magneton in inverse meters per tesla           2.542 623 58e-2       0.000 000 22e-2       m^-1 T^-1
nuclear magneton in K/T                                3.658 2637e-4         0.000 0064e-4         K T^-1
fine-structure constant                                7.297 352 568e-3      0.000 000 024e-3
inverse fine-structure constant                        137.035 999 11        0.000 000 46
Rydberg constant                                       10 973 731.568 525    0.000 073             m^-1
Rydberg constant times c in Hz                         3.289 841 960 360e15  0.000 000 000 022e15  Hz
Rydberg constant times hc in J                         2.179 872 09e-18      0.000 000 37e-18      J
Rydberg constant times hc in eV                        13.605 6923           0.000 0012            eV
Bohr radius                                            0.529 177 2108e-10    0.000 000 0018e-10    m
Hartree energy                                         4.359 744 17e-18      0.000 000 75e-18      J
Hartree energy in eV                                   27.211 3845           0.000 0023            eV
quantum of circulation                                 3.636 947 550e-4      0.000 000 024e-4      m^2 s^-1
quantum of circulation times 2                         7.273 895 101e-4      0.000 000 048e-4      m^2 s^-1
Fermi coupling constant                                1.166 39e-5           0.000 01e-5           GeV^-2
weak mixing angle                                      0.222 15              0.000 76
electron mass                                          9.109 3826e-31        0.000 0016e-31        kg
electron mass in u                                     5.485 799 0945e-4     0.000 000 0024e-4     u
electron mass energy equivalent                        8.187 1047e-14        0.000 0014e-14        J
electron mass energy equivalent in MeV                 0.510 998 918         0.000 000 044         MeV
electron-muon mass ratio                               4.836 331 67e-3       0.000 000 13e-3
electron-tau mass ratio                                2.875 64e-4           0.000 47e-4
electron-proton mass ratio                             5.446 170 2173e-4     0.000 000 0025e-4
electron-neutron mass ratio                            5.438 673 4481e-4     0.000 000 0038e-4
electron-deuteron mass ratio                           2.724 437 1095e-4     0.000 000 0013e-4
electron to alpha particle mass ratio                  1.370 933 555 75e-4   0.000 000 000 61e-4
electron charge to mass quotient                       -1.758 820 12e11      0.000 000 15e11       C kg^-1
electron molar mass                                    5.485 799 0945e-7     0.000 000 0024e-7     kg mol^-1
Compton wavelength                                     2.426 310 238e-12     0.000 000 016e-12     m
Compton wavelength over 2 pi                           386.159 2678e-15      0.000 0026e-15        m
classical electron radius                              2.817 940 325e-15     0.000 000 028e-15     m
Thomson cross section                                  0.665 245 873e-28     0.000 000 013e-28     m^2
electron magn. moment                                  -928.476 412e-26      0.000 080e-26         J T^-1
electron magn. moment to Bohr magneton ratio           -1.001 159 652 1859   0.000 000 000 0038
electron magn. moment to nuclear magneton ratio        -1838.281 971 07      0.000 000 85
electron magn. moment anomaly                          1.159 652 1859e-3     0.000 000 0038e-3
electron g factor                                      -2.002 319 304 3718   0.000 000 000 0075
electron-muon magn. moment ratio                       206.766 9894          0.000 0054
electron-proton magn. moment ratio                     -658.210 6862         0.000 0066
electron to shielded proton magn. moment ratio         -658.227 5956         0.000 0071
electron-neutron magn. moment ratio                    960.920 50            0.000 23
electron-deuteron magn. moment ratio                   -2143.923 493         0.000 023
electron to shielded helion magn. moment ratio         864.058 255           0.000 010
electron gyromagn. ratio                               1.760 859 74e11       0.000 000 15e11       s^-1 T^-1
electron gyromagn. ratio over 2 pi                     28 024.9532           0.0024                MHz T^-1
muon mass                                              1.883 531 40e-28      0.000 000 33e-28      kg
muon mass in u                                         0.113 428 9264        0.000 000 0030        u
muon mass energy equivalent                            1.692 833 60e-11      0.000 000 29e-11      J
muon mass energy equivalent in MeV                     105.658 3692          0.000 0094            MeV
muon-electron mass ratio                               206.768 2838          0.000 0054
muon-tau mass ratio                                    5.945 92e-2           0.000 97e-2
muon-proton mass ratio                                 0.112 609 5269        0.000 000 0029
muon-neutron mass ratio                                0.112 454 5175        0.000 000 0029
muon molar mass                                        0.113 428 9264e-3     0.000 000 0030e-3     kg mol^-1
muon Compton wavelength                                11.734 441 05e-15     0.000 000 30e-15      m
muon Compton wavelength over 2 pi                      1.867 594 298e-15     0.000 000 047e-15     m
muon magn. moment                                      -4.490 447 99e-26     0.000 000 40e-26      J T^-1
muon magn. moment to Bohr magneton ratio               -4.841 970 45e-3      0.000 000 13e-3
muon magn. moment to nuclear magneton ratio            -8.890 596 98         0.000 000 23
muon magn. moment anomaly                              1.165 919 81e-3       0.000 000 62e-3
muon g factor                                          -2.002 331 8396       0.000 000 0012
muon-proton magn. moment ratio                         -3.183 345 118        0.000 000 089
tau mass                                               3.167 77e-27          0.000 52e-27          kg
tau mass in u                                          1.907 68              0.000 31              u
tau mass energy equivalent                             2.847 05e-10          0.000 46e-10          J
tau mass energy equivalent in MeV                      1776.99               0.29                  MeV
tau-electron mass ratio                                3477.48               0.57
tau-muon mass ratio                                    16.8183               0.0027
tau-proton mass ratio                                  1.893 90              0.000 31
tau-neutron mass ratio                                 1.891 29              0.000 31
tau molar mass                                         1.907 68e-3           0.000 31e-3           kg mol^-1
tau Compton wavelength                                 0.697 72e-15          0.000 11e-15          m
tau Compton wavelength over 2 pi                       0.111 046e-15         0.000 018e-15         m
proton mass                                            1.672 621 71e-27      0.000 000 29e-27      kg
proton mass in u                                       1.007 276 466 88      0.000 000 000 13      u
proton mass energy equivalent                          1.503 277 43e-10      0.000 000 26e-10      J
proton mass energy equivalent in MeV                   938.272 029           0.000 080             MeV
proton-electron mass ratio                             1836.152 672 61       0.000 000 85
proton-muon mass ratio                                 8.880 243 33          0.000 000 23
proton-tau mass ratio                                  0.528 012             0.000 086
proton-neutron mass ratio                              0.998 623 478 72      0.000 000 000 58
proton charge to mass quotient                         9.578 833 76e7        0.000 000 82e7        C kg^-1
proton molar mass                                      1.007 276 466 88e-3   0.000 000 000 13e-3   kg mol^-1
proton Compton wavelength                              1.321 409 8555e-15    0.000 000 0088e-15    m
proton Compton wavelength over 2 pi                    0.210 308 9104e-15    0.000 000 0014e-15    m
proton magn. moment                                    1.410 606 71e-26      0.000 000 12e-26      J T^-1
proton magn. moment to Bohr magneton ratio             1.521 032 206e-3      0.000 000 015e-3
proton magn. moment to nuclear magneton ratio          2.792 847 351         0.000 000 028
proton g factor                                        5.585 694 701         0.000 000 056
proton-neutron magn. moment ratio                      -1.459 898 05         0.000 000 34
shielded proton magn. moment                           1.410 570 47e-26      0.000 000 12e-26      J T^-1
shielded proton magn. moment to Bohr magneton ratio    1.520 993 132e-3      0.000 000 016e-3
shielded proton magn. moment to nuclear magneton ratio 2.792 775 604         0.000 000 030
proton magn. shielding correction                      25.689e-6             0.015e-6
proton gyromagn. ratio                                 2.675 222 05e8        0.000 000 23e8        s^-1 T^-1
proton gyromagn. ratio over 2 pi                       42.577 4813           0.000 0037            MHz T^-1
shielded proton gyromagn. ratio                        2.675 153 33e8        0.000 000 23e8        s^-1 T^-1
shielded proton gyromagn. ratio over 2 pi              42.576 3875           0.000 0037            MHz T^-1
proton rms charge radius                               0.8750e-15            0.0068e-15            m
neutron mass                                           1.674 927 28e-27      0.000 000 29e-27      kg
neutron mass in u                                      1.008 664 915 60      0.000 000 000 55      u
neutron mass energy equivalent                         1.505 349 57e-10      0.000 000 26e-10      J
neutron mass energy equivalent in MeV                  939.565 360           0.000 081             MeV
neutron-electron mass ratio                            1838.683 6598         0.000 0013
neutron-muon mass ratio                                8.892 484 02          0.000 000 23
neutron-tau mass ratio                                 0.528 740             0.000 086
neutron-proton mass ratio                              1.001 378 418 70      0.000 000 000 58
neutron molar mass                                     1.008 664 915 60e-3   0.000 000 000 55e-3   kg mol^-1
neutron Compton wavelength                             1.319 590 9067e-15    0.000 000 0088e-15    m
neutron Compton wavelength over 2 pi                   0.210 019 4157e-15    0.000 000 0014e-15    m
neutron magn. moment                                   -0.966 236 45e-26     0.000 000 24e-26      J T^-1
neutron magn. moment to Bohr magneton ratio            -1.041 875 63e-3      0.000 000 25e-3
neutron magn. moment to nuclear magneton ratio         -1.913 042 73         0.000 000 45
neutron g factor                                       -3.826 085 46         0.000 000 90
neutron-electron magn. moment ratio                    1.040 668 82e-3       0.000 000 25e-3
neutron-proton magn. moment ratio                      -0.684 979 34         0.000 000 16
neutron to shielded proton magn. moment ratio          -0.684 996 94         0.000 000 16
neutron gyromagn. ratio                                1.832 471 83e8        0.000 000 46e8        s^-1 T^-1
neutron gyromagn. ratio over 2 pi                      29.164 6950           0.000 0073            MHz T^-1
deuteron mass                                          3.343 583 35e-27      0.000 000 57e-27      kg
deuteron mass in u                                     2.013 553 212 70      0.000 000 000 35      u
deuteron mass energy equivalent                        3.005 062 85e-10      0.000 000 51e-10      J
deuteron mass energy equivalent in MeV                 1875.612 82           0.000 16              MeV
deuteron-electron mass ratio                           3670.482 9652         0.000 0018
deuteron-proton mass ratio                             1.999 007 500 82      0.000 000 000 41
deuteron molar mass                                    2.013 553 212 70e-3   0.000 000 000 35e-3   kg mol^-1
deuteron magn. moment                                  0.433 073 482e-26     0.000 000 038e-26     J T^-1
deuteron magn. moment to Bohr magneton ratio           0.466 975 4567e-3     0.000 000 0050e-3
deuteron magn. moment to nuclear magneton ratio        0.857 438 2329        0.000 000 0092
deuteron-electron magn. moment ratio                   -4.664 345 548e-4     0.000 000 050e-4
deuteron-proton magn. moment ratio                     0.307 012 2084        0.000 000 0045
deuteron-neutron magn. moment ratio                    -0.448 206 52         0.000 000 11
deuteron rms charge radius                             2.1394e-15            0.0028e-15            m
helion mass                                            5.006 412 14e-27      0.000 000 86e-27      kg
helion mass in u                                       3.014 932 2434        0.000 000 0058        u
helion mass energy equivalent                          4.499 538 84e-10      0.000 000 77e-10      J
helion mass energy equivalent in MeV                   2808.391 42           0.000 24              MeV
helion-electron mass ratio                             5495.885 269          0.000 011
helion-proton mass ratio                               2.993 152 6671        0.000 000 0058
helion molar mass                                      3.014 932 2434e-3     0.000 000 0058e-3     kg mol^-1
shielded helion magn. moment                           -1.074 553 024e-26    0.000 000 093e-26     J T^-1
shielded helion magn. moment to Bohr magneton ratio    -1.158 671 474e-3     0.000 000 014e-3
shielded helion magn. moment to nuclear magneton ratio -2.127 497 723        0.000 000 025
shielded helion to proton magn. moment ratio           -0.761 766 562        0.000 000 012
shielded helion to shielded proton magn. moment ratio  -0.761 786 1313       0.000 000 0033
shielded helion gyromagn. ratio                        2.037 894 70e8        0.000 000 18e8        s^-1 T^-1
shielded helion gyromagn. ratio over 2 pi              32.434 1015           0.000 0028            MHz T^-1
alpha particle mass                                    6.644 6565e-27        0.000 0011e-27        kg
alpha particle mass in u                               4.001 506 179 149     0.000 000 000 056     u
alpha particle mass energy equivalent                  5.971 9194e-10        0.000 0010e-10        J
alpha particle mass energy equivalent in MeV           3727.379 17           0.000 32              MeV
alpha particle-electron mass ratio                     7294.299 5363         0.000 0032
alpha particle-proton mass ratio                       3.972 599 689 07      0.000 000 000 52
alpha particle molar mass                              4.001 506 179 149e-3  0.000 000 000 056e-3  kg mol^-1
Avogadro constant                                      6.022 1415e23         0.000 0010e23         mol^-1
atomic mass constant                                   1.660 538 86e-27      0.000 000 28e-27      kg
atomic mass constant energy equivalent                 1.492 417 90e-10      0.000 000 26e-10      J
atomic mass constant energy equivalent in MeV          931.494 043           0.000 080             MeV
Faraday constant                                       96 485.3383           0.0083                C mol^-1
Faraday constant for conventional electric current     96 485.336            0.016                 C_90 mol^-1
molar Planck constant                                  3.990 312 716e-10     0.000 000 027e-10     J s mol^-1
molar Planck constant times c                          0.119 626 565 72      0.000 000 000 80      J m mol^-1
molar gas constant                                     8.314 472             0.000 015             J mol^-1 K^-1
Boltzmann constant                                     1.380 6505e-23        0.000 0024e-23        J K^-1
Boltzmann constant in eV/K                             8.617 343e-5          0.000 015e-5          eV K^-1
Boltzmann constant in Hz/K                             2.083 6644e10         0.000 0036e10         Hz K^-1
Boltzmann constant in inverse meters per kelvin        69.503 56             0.000 12              m^-1 K^-1
molar volume of ideal gas (273.15 K, 101.325 kPa)      22.413 996e-3         0.000 039e-3          m^3 mol^-1
Loschmidt constant (273.15 K, 101.325 kPa)             2.686 7773e25         0.000 0047e25         m^-3
molar volume of ideal gas (273.15 K, 100 kPa)          22.710 981e-3         0.000 040e-3          m^3 mol^-1
Sackur-Tetrode constant (1 K, 100 kPa)                 -1.151 7047           0.000 0044
Sackur-Tetrode constant (1 K, 101.325 kPa)             -1.164 8677           0.000 0044
Stefan-Boltzmann constant                              5.670 400e-8          0.000 040e-8          W m^-2 K^-4
first radiation constant                               3.741 771 38e-16      0.000 000 64e-16      W m^2
first radiation constant for spectral radiance         1.191 042 82e-16      0.000 000 20e-16      W m^2 sr^-1
second radiation constant                              1.438 7752e-2         0.000 0025e-2         m K
Wien displacement law constant                         2.897 7685e-3         0.000 0051e-3         m K
molar mass of carbon-12                                12e-3                 (exact)               kg mol^-1
molar mass constant                                    1e-3                  (exact)               kg mol^-1
conventional value of Josephson constant               483 597.9e9           (exact)               Hz V^-1
conventional value of von Klitzing constant            25 812.807            (exact)               ohm
standard atmosphere                                    101 325               (exact)               Pa
standard acceleration of gravity                       9.806 65              (exact)               m s^-2
Cu x unit                                              1.002 077 10e-13      0.000 000 29e-13      m
Mo x unit                                              1.002 099 66e-13      0.000 000 53e-13      m
Angstrom star                                          1.000 015 09e-10      0.000 000 90e-10      m
lattice parameter of silicon                           543.102 122e-12       0.000 020e-12         m
{220} lattice spacing of silicon                       192.015 5965e-12      0.000 0070e-12        m
molar volume of silicon                                12.058 8382e-6        0.000 0024e-6         m^3 mol^-1
electron volt                                          1.602 176 53e-19      0.000 000 14e-19      J
unified atomic mass unit                               1.660 538 86e-27      0.000 000 28e-27      kg
natural unit of velocity                               299 792 458           (exact)               m s^-1
natural unit of action                                 1.054 571 68e-34      0.000 000 18e-34      J s
natural unit of action in eV s                         6.582 119 15e-16      0.000 000 56e-16      eV s
natural unit of mass                                   9.109 3826e-31        0.000 0016e-31        kg
natural unit of energy                                 8.187 1047e-14        0.000 0014e-14        J
natural unit of energy in MeV                          0.510 998 918         0.000 000 044         MeV
natural unit of momentum                               2.730 924 19e-22      0.000 000 47e-22      kg m s^-1
natural unit of momentum in MeV/c                      0.510 998 918         0.000 000 044         MeV/c
natural unit of length                                 386.159 2678e-15      0.000 0026e-15        m
natural unit of time                                   1.288 088 6677e-21    0.000 000 0086e-21    s
atomic unit of charge                                  1.602 176 53e-19      0.000 000 14e-19      C
atomic unit of mass                                    9.109 3826e-31        0.000 0016e-31        kg
atomic unit of action                                  1.054 571 68e-34      0.000 000 18e-34      J s
atomic unit of length                                  0.529 177 2108e-10    0.000 000 0018e-10    m
atomic unit of energy                                  4.359 744 17e-18      0.000 000 75e-18      J
atomic unit of time                                    2.418 884 326 505e-17 0.000 000 000 016e-17 s
atomic unit of force                                   8.238 7225e-8         0.000 0014e-8         N
atomic unit of velocity                                2.187 691 2633e6      0.000 000 0073e6      m s^-1
atomic unit of momentum                                1.992 851 66e-24      0.000 000 34e-24      kg m s^-1
atomic unit of current                                 6.623 617 82e-3       0.000 000 57e-3       A
atomic unit of charge density                          1.081 202 317e12      0.000 000 093e12      C m^-3
atomic unit of electric potential                      27.211 3845           0.000 0023            V
atomic unit of electric field                          5.142 206 42e11       0.000 000 44e11       V m^-1
atomic unit of electric field gradient                 9.717 361 82e21       0.000 000 83e21       V m^-2
atomic unit of electric dipole moment                  8.478 353 09e-30      0.000 000 73e-30      C m
atomic unit of electric quadrupole moment              4.486 551 24e-40      0.000 000 39e-40      C m^2
atomic unit of electric polarizablity                  1.648 777 274e-41     0.000 000 016e-41     C^2 m^2 J^-1
atomic unit of 1st hyperpolarizablity                  3.206 361 51e-53      0.000 000 28e-53      C^3 m^3 J^-2
atomic unit of 2nd hyperpolarizablity                  6.235 3808e-65        0.000 0011e-65        C^4 m^4 J^-3
atomic unit of magn. flux density                      2.350 517 42e5        0.000 000 20e5        T
atomic unit of magn. dipole moment                     1.854 801 90e-23      0.000 000 16e-23      J T^-1
atomic unit of magnetizability                         7.891 036 60e-29      0.000 000 13e-29      J T^-2
atomic unit of permittivity                            1.112 650 056...e-10  (exact)               F m^-1
joule-kilogram relationship                            1.112 650 056...e-17  (exact)               kg
joule-inverse meter relationship                       5.034 117 20e24       0.000 000 86e24       m^-1
joule-hertz relationship                               1.509 190 37e33       0.000 000 26e33       Hz
joule-kelvin relationship                              7.242 963e22          0.000 013e22          K
joule-electron volt relationship                       6.241 509 47e18       0.000 000 53e18       eV
joule-atomic mass unit relationship                    6.700 5361e9          0.000 0011e9          u
joule-hartree relationship                             2.293 712 57e17       0.000 000 39e17       E_h
kilogram-joule relationship                            8.987 551 787...e16   (exact)               J
kilogram-inverse meter relationship                    4.524 438 91e41       0.000 000 77e41       m^-1
kilogram-hertz relationship                            1.356 392 66e50       0.000 000 23e50       Hz
kilogram-kelvin relationship                           6.509 650e39          0.000 011e39          K
kilogram-electron volt relationship                    5.609 588 96e35       0.000 000 48e35       eV
kilogram-atomic mass unit relationship                 6.022 1415e26         0.000 0010e26         u
kilogram-hartree relationship                          2.061 486 05e34       0.000 000 35e34       E_h
inverse meter-joule relationship                       1.986 445 61e-25      0.000 000 34e-25      J
inverse meter-kilogram relationship                    2.210 218 81e-42      0.000 000 38e-42      kg
inverse meter-hertz relationship                       299 792 458           (exact)               Hz
inverse meter-kelvin relationship                      1.438 7752e-2         0.000 0025e-2         K
inverse meter-electron volt relationship               1.239 841 91e-6       0.000 000 11e-6       eV
inverse meter-atomic mass unit relationship            1.331 025 0506e-15    0.000 000 0089e-15    u
inverse meter-hartree relationship                     4.556 335 252 760e-8  0.000 000 000 030e-8  E_h
hertz-joule relationship                               6.626 0693e-34        0.000 0011e-34        J
hertz-kilogram relationship                            7.372 4964e-51        0.000 0013e-51        kg
hertz-inverse meter relationship                       3.335 640 951...e-9   (exact)               m^-1
hertz-kelvin relationship                              4.799 2374e-11        0.000 0084e-11        K
hertz-electron volt relationship                       4.135 667 43e-15      0.000 000 35e-15      eV
hertz-atomic mass unit relationship                    4.439 821 667e-24     0.000 000 030e-24     u
hertz-hartree relationship                             1.519 829 846 006e-16 0.000 000 000 010e-16 E_h
kelvin-joule relationship                              1.380 6505e-23        0.000 0024e-23        J
kelvin-kilogram relationship                           1.536 1808e-40        0.000 0027e-40        kg
kelvin-inverse meter relationship                      69.503 56             0.000 12              m^-1
kelvin-hertz relationship                              2.083 6644e10         0.000 0036e10         Hz
kelvin-electron volt relationship                      8.617 343e-5          0.000 015e-5          eV
kelvin-atomic mass unit relationship                   9.251 098e-14         0.000 016e-14         u
kelvin-hartree relationship                            3.166 8153e-6         0.000 0055e-6         E_h
electron volt-joule relationship                       1.602 176 53e-19      0.000 000 14e-19      J
electron volt-kilogram relationship                    1.782 661 81e-36      0.000 000 15e-36      kg
electron volt-inverse meter relationship               8.065 544 45e5        0.000 000 69e5        m^-1
electron volt-hertz relationship                       2.417 989 40e14       0.000 000 21e14       Hz
electron volt-kelvin relationship                      1.160 4505e4          0.000 0020e4          K
electron volt-atomic mass unit relationship            1.073 544 171e-9      0.000 000 092e-9      u
electron volt-hartree relationship                     3.674 932 45e-2       0.000 000 31e-2       E_h
atomic mass unit-joule relationship                    1.492 417 90e-10      0.000 000 26e-10      J
atomic mass unit-kilogram relationship                 1.660 538 86e-27      0.000 000 28e-27      kg
atomic mass unit-inverse meter relationship            7.513 006 608e14      0.000 000 050e14      m^-1
atomic mass unit-hertz relationship                    2.252 342 718e23      0.000 000 015e23      Hz
atomic mass unit-kelvin relationship                   1.080 9527e13         0.000 0019e13         K
atomic mass unit-electron volt relationship            931.494 043e6         0.000 080e6           eV
atomic mass unit-hartree relationship                  3.423 177 686e7       0.000 000 023e7       E_h
hartree-joule relationship                             4.359 744 17e-18      0.000 000 75e-18      J
hartree-kilogram relationship                          4.850 869 60e-35      0.000 000 83e-35      kg
hartree-inverse meter relationship                     2.194 746 313 705e7   0.000 000 000 015e7   m^-1
hartree-hertz relationship                             6.579 683 920 721e15  0.000 000 000 044e15  Hz
hartree-kelvin relationship                            3.157 7465e5          0.000 0055e5          K
hartree-electron volt relationship                     27.211 3845           0.000 0023            eV
hartree-atomic mass unit relationship                  2.921 262 323e-8      0.000 000 019e-8      u

*/

constexpr double PC_CGS_L         = 1.0e+02;                           // CGS length    per SI unit
constexpr double PC_CGS_M         = 1.0e+03;                           // CGS mass      per SI unit
constexpr double PC_CGS_T         = 1.0e+00;                           // CGS time      per SI unit

constexpr double PC_CGS_L2        = (PC_CGS_L*PC_CGS_L);               // CGS length^2  per SI unit
constexpr double PC_CGS_L3        = (PC_CGS_L2*PC_CGS_L);              // CGS length^3  per SI unit
constexpr double PC_CGS_M2        = (PC_CGS_M*PC_CGS_M);               // CGS mass^2    per SI unit
constexpr double PC_CGS_M3        = (PC_CGS_M2*PC_CGS_M);              // CGS mass^3    per SI unit
constexpr double PC_CGS_T2        = (PC_CGS_T*PC_CGS_T);               // CGS time^2    per SI unit
constexpr double PC_CGS_T3        = (PC_CGS_T2*PC_CGS_T);              // CGS time^3    per SI unit
constexpr double PC_CGS_V         = (PC_CGS_L/PC_CGS_T);               // CGS velocity  per SI unit
constexpr double PC_CGS_F         = (PC_CGS_M*PC_CGS_L/PC_CGS_T2);     // CGS force     per SI unit
constexpr double PC_CGS_E         = (PC_CGS_F*PC_CGS_L);               // CGS energy    per SI unit
constexpr double PC_CGS_Q         = 2.99792458e+09;                    // CGS charge    per SI unit - |c/(m/s)|*10
constexpr double PC_CGS_FLUX      = (1.0/PC_CGS_L2/PC_CGS_T);          // CGS flux      per SI unit - (i.e. 1/cm^2/s)
constexpr double PC_CGS_E_FLUX    = (PC_CGS_E/PC_CGS_L2/PC_CGS_T);     // CGS energy fl per SI unit - (i.e. erg/cm^2/s)

constexpr double PC_SI_C          = 2.99792458e+08;                    // vacuum speed of light              [m/s]
constexpr double PC_SI_MU0        = (4.0*M_PI*1e-7);                   // magnetic flux                    [N/A^2]
constexpr double PC_SI_EPSILON0   = (1.0/(PC_SI_MU0*PC_SI_C*PC_SI_C)); // electric flux                      [F/m]
constexpr double PC_SI_C_ERR      = 0.00000000e+08;                    // vacuum speed of light              [m/s]
constexpr double PC_SI_H          = 6.6260693e-34;                     // Plank                              [J s]
constexpr double PC_SI_H_ERR      = 0.0000011e-34;                     // Plank                              [J s]
constexpr double PC_SI_HBAR       = (PC_SI_H/2.0/M_PI);                // h / 2pi                            [J s]
constexpr double PC_SI_HBAR_ERR   = (PC_SI_H_ERR/2.0/M_PI);            // h / 2pi                            [J s]
constexpr double PC_SI_G          = 6.6742e-11;                        // grav                        [m^3/kg/s^2]
constexpr double PC_SI_G_ERR      = 0.0010e-11;                        // grav                        [m^3/kg/s^2]
constexpr double PC_SI_E          = 1.60217653e-19;                    // electric charge unit                 [C]
constexpr double PC_SI_E_ERR      = 0.00000014e-19;                    // electric charge unit                 [C]
constexpr double PC_SI_K          = 1.3806505e-23;                     // Boltzmann                          [J/K]
constexpr double PC_SI_K_ERR      = 0.0000024e-23;                     // Boltzmann                          [J/K]
constexpr double PC_SI_SIGMA      = 5.670400e-8;                       // Stefan-Boltzmann             [W/m^2/K^4]
constexpr double PC_SI_SIGMA_ERR  = 0.000040e-8;                       // pi^2 k^4/60 hbar^3 c^2       [W/m^2/K^4]

constexpr double PC_SI_AMU        = 1.66053886e-27;                    // amu-mass                            [kg]
constexpr double PC_SI_AMU_ERR    = 0.00000028e-27;                    // amu-mass                            [kg]
constexpr double PC_SI_EV         = 1.60217653e-19;                    // 1 eV                                 [J]
constexpr double PC_SI_EV_ERR     = 0.00000014e-19;                    // 1 eV                                 [J]

constexpr double PC_SI_E_MASS     =  9.1093826e-31;                    // e-mass                              [kg]
constexpr double PC_SI_E_MASS_ERR = 0.0000016e-28;                     // e-mass                              [kg]
constexpr double PC_SI_P_MASS     = 1.67262171e-27;                    // p-mass                              [kg]
constexpr double PC_SI_P_MASS_ERR = 0.00000029e-24;                    // p-mass                              [kg]

constexpr double PC_SI_TROPICAL_YEAR = 3.15569260e+07;                 // tropical year                        [s]
constexpr double PC_SI_JULIAN_YEAR   = 3.15576e+07;                    // julian year                          [s]

constexpr double PC_SI_AU            = 1.4959787e+11;                  // astronomical unit                    [m]
constexpr double PC_SI_LIGHTYEAR     = (PC_SI_C*PC_SI_JULIAN_YEAR);    // light year                           [m]
constexpr double PC_SI_PARSEC        = 3.0856776e+16;                  // parsec                               [m]
constexpr double PC_SI_SUN_MASS      = 1.9891e+30;                     // solar mass                          [kg]
constexpr double PC_SI_SUN_RADIUS    = 6.96e+8;                        // solar radius                         [m]
constexpr double PC_SI_SUN_LUMINOSITY= 3.9e+26;                        // solar luminosity                     [W]

constexpr double PC_G4_L          = 1.0e+02;                           // G4 length     per SI unit
constexpr double PC_G4_T          = 1.0e+09;                           // G4 time       per SI unit
constexpr double PC_G4_V          = (PC_G4_L/PC_G4_T);                 // G4 velocity   per SI unit

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Some useful numerical constants
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
constexpr double num_pi             = M_PI;                        // pi
constexpr double num_e              = M_E;                         // base of natural log
constexpr double num_ln_10          = M_LN10;                      // natural log of 10
constexpr double num_log10_e        = M_LOG10E;                    // log base 10 of e
constexpr double num_root_two       = M_SQRT2;                     // sqrt(2)
constexpr double num_root_half      = M_SQRT1_2;                   // sqrt(0.5)

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Unit conversion
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
constexpr double cgs_L              = PC_CGS_L;                    // CGS length    per SI unit
constexpr double cgs_M              = PC_CGS_M;                    // CGS mass      per SI unit
constexpr double cgs_T              = PC_CGS_T;                    // CGS time      per SI unit

constexpr double cgs_L2             = PC_CGS_L2;                   // CGS length^2  per SI unit
constexpr double cgs_L3             = PC_CGS_L3;                   // CGS length^3  per SI unit
constexpr double cgs_M2             = PC_CGS_M2;                   // CGS mass^2    per SI unit
constexpr double cgs_M3             = PC_CGS_M3;                   // CGS mass^3    per SI unit
constexpr double cgs_T2             = PC_CGS_T2;                   // CGS time^2    per SI unit
constexpr double cgs_T3             = PC_CGS_T3;                   // CGS time^3    per SI unit
constexpr double cgs_V              = PC_CGS_V;                    // CGS velocity  per SI unit
constexpr double cgs_F              = PC_CGS_F;                    // CGS force     per SI unit
constexpr double cgs_E              = PC_CGS_E;                    // CGS energy    per SI unit
constexpr double cgs_Q              = PC_CGS_Q;                    // CGS charge    per SI unit - |c/(m/s)|*10
constexpr double cgs_FLUX           = PC_CGS_FLUX;                 // CGS flux      per SI unit - (i.e. 1/cm^2/s)
constexpr double cgs_E_FLUX         = PC_CGS_E_FLUX;               // CGS energy fl per SI unit - (i.e. erg/cm^2/s)

constexpr double g4_L               = PC_G4_L;                     // G4 length     per SI unit
constexpr double g4_T               = PC_G4_T;                     // G4 time       per SI unit
constexpr double g4_V               = PC_G4_V;                     // G4 velocity   per SI unit

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constants and errors from NIST, Mohr and Taylor, 2002 (as included in comments above)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
constexpr double num_avogadro       = 6.0221415e+23;               // avogadro                              [1/mol]
constexpr double num_avogadro_err   = 0.0000010e+23;               // avogadro                              [1/mol]
constexpr double num_alpha          = 7.297352568e-3;              // fine structure                            [1]
constexpr double num_alpha_err      = 0.000000024e-3;              // fine structure                            [1]

constexpr double si_c               = PC_SI_C;                     // vacuum speed of light                   [m/s]
constexpr double si_1_c             = 1.0/si_c;                    // reciprocal of vacuum speed of light     [s/m]
constexpr double si_mu0             = PC_SI_MU0;                   // magnetic flux                         [N/A^2]
constexpr double si_epsilon0        = PC_SI_EPSILON0;              // electric flux                           [F/m]
constexpr double si_c_err           = PC_SI_C_ERR;                 // vacuum speed of light                   [m/s]
constexpr double si_h               = PC_SI_H;                     // Plank                                   [J s]
constexpr double si_h_err           = PC_SI_H_ERR;                 // Plank                                   [J s]
constexpr double si_hbar            = PC_SI_HBAR;                  // h / 2pi                                 [J s]
constexpr double si_hbar_err        = PC_SI_HBAR_ERR;              // h / 2pi                                 [J s]
constexpr double si_G               = PC_SI_G;                     // grav                             [m^3/kg/s^2]
constexpr double si_G_err           = PC_SI_G_ERR;                 // grav                             [m^3/kg/s^2]
constexpr double si_e               = PC_SI_E;                     // electric charge unit                      [C]
constexpr double si_e_err           = PC_SI_E_ERR;                 // electric charge unit                      [C]
constexpr double si_k               = PC_SI_K;                     // Boltzmann                               [J/K]
constexpr double si_k_err           = PC_SI_K_ERR;                 // Boltzmann                               [J/K]
constexpr double si_sigma           = PC_SI_SIGMA;                 // Stefan-Boltzmann                  [W/m^2/K^4]
constexpr double si_sigma_err       = PC_SI_SIGMA_ERR;             // pi^2 k^4/60 hbar^3 c^2            [W/m^2/K^4]

constexpr double cgs_c              = PC_SI_C*PC_CGS_V;            // vacuum speed of light                  [cm/s]
constexpr double cgs_c_err          = PC_SI_C_ERR*PC_CGS_V;        // vacuum speed of light                  [cm/s]
constexpr double cgs_1_c            = 1.0/cgs_c;                   // reciprocal of vacuum speed of light    [s/cm]
constexpr double cgs_h              = PC_SI_H*PC_CGS_E/PC_CGS_T;   // Plank                                 [erg s]
constexpr double cgs_h_err          = PC_SI_H_ERR*PC_CGS_E/PC_CGS_T;// Plank                                [erg s]
constexpr double cgs_hbar           = PC_SI_H*PC_CGS_E/PC_CGS_T;   // h / 2pi                               [erg s]
constexpr double cgs_hbar_err       = PC_SI_H_ERR*PC_CGS_E/PC_CGS_T;// h / 2pi                              [erg s]
constexpr double cgs_G              = PC_SI_G*PC_CGS_L3/PC_CGS_M/PC_CGS_T2; // grav                    [cm^3/g/s^2]
constexpr double cgs_G_err          = PC_SI_G_ERR*PC_CGS_L3/PC_CGS_M/PC_CGS_T2; //  grav               [cm^3/g/s^2]
constexpr double cgs_e              = PC_SI_E*PC_CGS_Q;            // electric charge unit                    [esu]
constexpr double cgs_e_err          = PC_SI_E_ERR*PC_CGS_Q;        // electric charge unit                    [esu]
constexpr double cgs_k              = PC_SI_K*PC_CGS_E;            // Boltzmann                             [erg/K]
constexpr double cgs_k_err          = PC_SI_K_ERR*PC_CGS_E;        // Boltzmann                             [erg/K]
constexpr double cgs_sigma          = PC_SI_SIGMA*PC_CGS_E/PC_CGS_L2; // Stefan-Boltzmann            [erg/cm^2/K^4]
constexpr double cgs_sigma_err      = PC_SI_SIGMA_ERR*PC_CGS_E/PC_CGS_L2; // pi^2 k^4/60 hbar^3 c^2      [as above]

constexpr double g4_c               = PC_SI_C*PC_G4_V;             // vacuum speed of light                  [cm/s]
constexpr double g4_1_c             = 1.0/g4_c;                    // reciprocal of vacuum speed of light    [s/cm]

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Unit conversion constants
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
constexpr double si_amu             = PC_SI_AMU;                   // amu-mass                                 [kg]
constexpr double si_amu_err         = PC_SI_AMU_ERR;               // amu-mass                                 [kg]
constexpr double si_eV              = PC_SI_EV;                    // 1 eV                                      [J]
constexpr double si_eV_err          = PC_SI_EV_ERR;                // 1 eV                                      [J]

constexpr double cgs_amu            = PC_SI_AMU*PC_CGS_M;          // amu-mass                                  [g]
constexpr double cgs_amu_err        = PC_SI_AMU_ERR*PC_CGS_M;      // amu-mass                                  [g]
constexpr double cgs_eV             = PC_SI_EV*PC_CGS_E;           // 1 eV                                    [erg]
constexpr double cgs_eV_err         = PC_SI_EV_ERR*PC_CGS_E;       // 1 eV                                    [erg]

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Particle properties
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
constexpr double si_e_chg           = -PC_SI_E;                    // e-charge                                  [C]
constexpr double si_e_chg_err       = PC_SI_E_ERR;                 // e-charge                                  [C]
constexpr double si_e_mass          = PC_SI_E_MASS;                // e-mass                                   [kg]
constexpr double si_e_mass_err      = PC_SI_E_MASS_ERR;            // e-mass                                   [kg]
constexpr double si_p_chg           = PC_SI_E;                     // e-charge                                  [C]
constexpr double si_p_chg_err       = PC_SI_E_ERR;                 // e-charge                                  [C]
constexpr double si_p_mass          = PC_SI_P_MASS;                // p-mass                                   [kg]
constexpr double si_p_mass_err      = PC_SI_P_MASS_ERR;            // p-mass                                   [kg]

constexpr double cgs_e_chg          = -PC_SI_E*PC_CGS_Q;           // e-charge                                [esu]
constexpr double cgs_e_chg_err      = PC_SI_E_ERR*PC_CGS_Q;        // e-charge                                [esu]
constexpr double cgs_e_mass         = PC_SI_E_MASS*PC_CGS_M;       // e-mass                                    [g]
constexpr double cgs_e_mass_err     = PC_SI_E_MASS_ERR*PC_CGS_M;   // e-mass                                    [g]
constexpr double cgs_p_chg          = PC_SI_E*PC_CGS_Q;            // e-charge                                [esu]
constexpr double cgs_p_chg_err      = PC_SI_E_ERR*PC_CGS_Q;        // e-charge                                [esu]
constexpr double cgs_p_mass         = PC_SI_P_MASS*PC_CGS_M;       // p-mass                                    [g]
constexpr double cgs_p_mass_err     = PC_SI_P_MASS_ERR*PC_CGS_M;   // p-mass                                    [g]

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Astronomical constants
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
constexpr double si_tropical_year   = PC_SI_TROPICAL_YEAR;         // tropical year                             [s]
constexpr double si_julian_year     = PC_SI_JULIAN_YEAR;           // julian year                               [s]

constexpr double cgs_tropical_year  = PC_SI_TROPICAL_YEAR*PC_CGS_T;// tropical year                             [s]
constexpr double cgs_julian_year    = PC_SI_JULIAN_YEAR*PC_CGS_T;  // julian year                               [s]

constexpr double si_au              = PC_SI_AU;                    // astronomical unit                         [m]
constexpr double si_light_year      = PC_SI_LIGHTYEAR;             // light year                                [m]
constexpr double si_parsec          = PC_SI_PARSEC;                // parsec                                    [m]
constexpr double si_sun_mass        = PC_SI_SUN_MASS;              // solar mass                               [kg]
constexpr double si_sun_radius      = PC_SI_SUN_RADIUS;            // solar radius                              [m]
constexpr double si_sun_luminosity  = PC_SI_SUN_LUMINOSITY;        // solar luminosity                          [W]

constexpr double cgs_au             = PC_SI_AU*PC_CGS_L;           // astrnomical unit                         [cm]
constexpr double cgs_light_year     = PC_SI_LIGHTYEAR*PC_CGS_L;    // light year                               [cm]
constexpr double cgs_parsec         = PC_SI_PARSEC*PC_CGS_L;       // parsec                                   [cm]
constexpr double cgs_sun_mass       = PC_SI_SUN_MASS*PC_CGS_M;     // solar mass                                [g]
constexpr double cgs_sun_radius     = PC_SI_SUN_RADIUS*PC_CGS_L;   // solar radius                             [cm]
constexpr double cgs_sun_luminosity = PC_SI_SUN_LUMINOSITY*PC_CGS_E/PC_CGS_T; // solar luminosity           [erg/s]

#if 0
template<class T> std::string convertToString(const T& x)
{
  std::ostringstream stream;
  stream << x;
  return stream.str();
}

template<> std::string convertToString<>(const std::string& x);
template<> std::string convertToString<>(const float& x);
template<> std::string convertToString<>(const double& x);
template<> std::string convertToString<>(const long double& x);

template<class T> bool convertFromString(const std::string& str, T& x)
{
  std::istringstream stream(str);
  return(stream >> x);
}

template<> bool convertFromString<>(const std::string& str, std::string& x);

template<class T> std::string
FDAV(const std::string desc, const T& val, const std::string units="",
     unsigned min_desc_len=0, unsigned indent=0)
{
  std::ostringstream s;
  s << std::string(indent*2,' ');
  if(min_desc_len != 0)s << std::setw(min_desc_len-indent*2) << std::left;
  s << desc << ' ' << convertToString(val);
  if(units.length() != 0)s << ' ' << units;
  return s.str();
}

template<class T> void tokenize(const std::string& csl, std::set<T>& list)
{
  static const char sep[] = ", \t\r\n";
  std::string::size_type next = 0;
  while((next != std::string::npos)&&(next < csl.length()))
  {
    std::string::size_type here = next;
    next = csl.find_first_of(sep,here);
    std::string val_string;
    if(next == std::string::npos)val_string = csl.substr(here);
    else val_string = csl.substr(here,(next++)-here);
    if(!val_string.empty())
    {
      std::istringstream val_stream(val_string);
      T val;
      if(val_stream >> val)list.insert(val);
    }
  }
}
#endif

} } } // namespace calin::math::constants
