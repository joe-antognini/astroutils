#! /usr/bin/env python

#
# Cosmology utilities
#

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
    return (1 / sqrt(omega_k - 1) * sin(i_integral(redshift, omega_m,
      omega_l, omega_r, omega_k) * sqrt(omega_k - 1)))
  elif omega_k > 1:
    return (1 / sqrt(omega_k - 1) * sinh(i_integral(redshift, omega_m,
      omega_l, omega_r, omega_k) * sqrt(omega_k - 1)))

def comove_coord(z, omega_m=.27, omega_l=.73, omega_r=0, omega_k=1):
  '''Calculates the comoving coordinate of eq. 29.176 of Carroll & Ostlie
  pg. 1210.  Returns the coordinate in Mpc.'''
  import consts
  return (consts.speed_of_light_km_s / consts.hubble_constant * 
    s_func(z, omega_m, omega_l, omega_r, omega_k))

def luminosity_distance(z, omega_m=.27, omega_l=.73, omega_r=0, omega_k=1):
  '''Returns the luminosity distance of a particular redshift.  One can
  optionally change the cosmology.  Otherwise it assumes a standard
  concordance cosmology.  The distance is returned in Mpc.'''

  return comove_coord(z, omega_m, omega_l, omega_r, omega_k) * (1 + z)

def angular_distance(z, omega_m=.27, omega_l=.73, omega_r=0, omega_k=1):
  '''Returns the angular diameter distance of a particular redshift.  One
  can optionally change the cosmology.  Otherwise it assumes a standard
  concordance cosmology.  The distance is returned in Mpc.'''

  return (luminosity_distance(z, omega_m, omega_l, omega_r, omega_k) /
    (1 + z)**2)
   
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

