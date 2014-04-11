#! /usr/bin/env python

#
# Dynamics routines
#

from numpy import sqrt
from astropy.constants import G, M_sun

def vcrit_tripsing(m, a):
  '''Calculate the critical velocity of a single star scattering off of a
  hierarchical triple.  (The critical velocity is the velocity at which the
  total energy of the system is zero.)

  Inputs:
    m -- a tuple containing the four masses of the stars as astropy
      quantities in the following order:
        (star 1 of the inner binary, star 2 of the inner binary, tertiary,
        interloping star)

    a -- a tuple containing the two semi-major axes of the hierarchical
      triple as an astropy quantity in the following order:
      (inner semi-major axis, outer semi-major axis)

    Output:
      The critical velocity as an astropy quantity.
  '''

  # Unpack the inputs
  m1, m2, m3, m4 = m
  a11, a1 = a

  return sqrt((G * (a1 * m1 * m2 + a11 * (m1 + m2) * m3) * (m1 + m2 + m3 +
    m4)) / (a1 * a11 * (m1 + m2 + m3) * m4))

def mardling(q_out, e_out, inc=0):
  '''Calculate the maximum eccentricity for which the triple is stable
  according to the Mardling stability criterion.  See Mardling & Aarseth
  (2001) for more details, specifically Eq. 90.

  Inputs:
    q_out -- The mass ratio between the inner and outer binaries.
      Equivalent to m3 / (m1 + m2)

    e_out -- Eccentricity of the outer orbit.

    inc [deg ] -- Inclination of the outer binary.  

  Output:
    The maximum eccentricity for which the triple is stable.
  '''
  from math import sqrt
  
  # An empirically determined constant
  C = 2.8

  return 1 - C * ((1 + q_out) * (1 + e_out) / sqrt(1 - e_out))**(2./5)
