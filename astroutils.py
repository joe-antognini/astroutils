#! /usr/bin/env python -W ignore::DeprecationWarning

#
# astroutils
#
# A collection of functions useful for various astronomical purposes.
#
# Joe Antognini
# Sun Mar 13 20:16:21 EDT 2011
#

###
### CONSTANTS
###

hubble_constant = 70.4 # km / s / Mpc
speed_of_light_km_s = 299792.458 # km / s
solarmass2kg = 1.989e30
parsec2m = 3.0857e16
au2m = 1.496e11
year2sec = 3.155693e7
rsun2m = 6.9599e8

def mpc2meter(dist):
  '''Converts Megaparsecs to meters.'''
  return 3.0857e+22 * dist

def kpc2meter(dist):
  '''Converts kiloparsecs to meters.'''
  return 3.0857e+19 * dist

def meter2mpc(dist):
  '''Converts meters to Megaparsecs.'''
  return dist / 3.0857e+22

def meter2kpc(dist):
  '''Converts meters to kiloparsecs.'''
  return dist / 3.0857e+19

def radian2arcsec(ang):
  '''Converts an angle in radians to arcseconds.'''
  return ang * 206265.

def arcsec2radian(ang):
  '''Converts an angle in arcseconds to radians.'''
  return ang / 206265.

def radian2arcmin(ang):
  '''Converts an angle in radians to arcminutes.'''
  return ang * 3437.75

def arcmin2radian(ang):
  '''Converts an angle in arcminutes to radians.'''
  return ang / 3437.75

def i_integral(redshift, omega_m=.27, omega_l=.73, omega_r=0, omega_k=1):
  '''Calculates the integral of Carroll & Ostlie pg. 1209.'''
  from scipy.integrate import quad
  from math import sqrt

  return quad(lambda z: 1 / sqrt(omega_m * (1 + z)**3 + omega_r * \
    (1 + z)**4 + omega_l + (1 - omega_k) * (1 + z)**2), 0, redshift)[0]

def s_func(redshift, omega_m=.27, omega_l=.73, omega_r=0, omega_k=1):
  '''Calculates 29.173-175 of Carroll & Ostlie pg. 1210'''
  from math import sin, sinh, sqrt
  if omega_k == 1:
    return i_integral(redshift, omega_m, omega_l, omega_r, omega_k)
  elif omega_k < 1:
    return 1 / sqrt(omega_k - 1) * sin(i_integral(redshift, omega_m, \
      omega_l, omega_r, omega_k) * sqrt(omega_k - 1))
  elif omega_k > 1:
    return 1 / sqrt(omega_k - 1) * sinh(i_integral(redshift, omega_m, \
      omega_l, omega_r, omega_k) * sqrt(omega_k - 1))

def comove_coord(z, omega_m=.27, omega_l=.73, omega_r=0, omega_k=1):
  '''Calculates the comoving coordinate of eq. 29.176 of Carroll & Ostlie
  pg. 1210.  Returns the coordinate in Mpc.'''
  return speed_of_light_km_s / hubble_constant * s_func(z, omega_m, \
    omega_l, omega_r, omega_k)

def luminosity_distance(z, omega_m=.27, omega_l=.73, omega_r=0, omega_k=1):
  '''Returns the luminosity distance of a particular redshift.  One can
  optionally change the cosmology.  Otherwise it assumes a standard
  concordance cosmology.  The distance is returned in Mpc.'''

  return comove_coord(z, omega_m, omega_l, omega_r, omega_k) * (1 + z)

def angular_distance(z, omega_m=.27, omega_l=.73, omega_r=0, omega_k=1):
  '''Returns the angular diameter distance of a particular redshift.  One
  can optionally change the cosmology.  Otherwise it assumes a standard
  concordance cosmology.  The distance is returned in Mpc.'''

  return luminosity_distance(z, omega_m, omega_l, omega_r, omega_k) / \
    (1 + z)**2
   
def universe_age(redshift, omega_m=.27, omega_l=.73, omega_r=8.24e-5, \
  omega_k=1):
  '''Returns the age of the universe at a given redshift.  Assumes a
  standard lambda CDM cosmology, but the cosmology can be optionally
  changed.  The age is returned in Myr.'''

  from scipy.integrate import quad, inf
  from math import sqrt

  h_const = hubble_constant * 1.022e-6 # Convert H_0 to 1/Myr
  return 1 / h_const * quad(lambda z: 1 / ((1 + z) * sqrt(omega_m * \
    (1 + z)**3 + omega_l + omega_r * (1 + z)**4 + (1 - omega_k) * \
    (1 + z)**2)), redshift, inf)[0]

def manom_eanom(e, m):
  '''Calculate the eccentric anomaly given the eccentricity and mean
  anomaly of an orbit.'''
  from scipy.optimize import brentq
  from math import sin, pi
  
  if m < 0 or m > 2 * pi:
    raise ValueError('mean anomaly must be an angle between 0 and 2pi')

  return brentq(lambda x: x - e * sin(x) - m, 0, 2 * pi)
