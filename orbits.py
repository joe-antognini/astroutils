#! /usr/bin/env python

#
# Orbital mechanics utilities
#

def manom_eanom(e, m):
  '''Calculate the eccentric anomaly given the eccentricity and mean
  anomaly of an orbit.'''
  from scipy.optimize import brentq
  from math import sin, pi
  
  if m < 0 or m > 2 * pi:
    raise ValueError('mean anomaly must be an angle between 0 and 2pi')

  return brentq(lambda x: x - e * sin(x) - m, 0, 2 * pi)

