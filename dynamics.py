#! /usr/bin/env python

#
# Dynamics routines
#

from numpy import sqrt

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

  from astropy.constants import G

  # Unpack the inputs
  m1, m2, m3, m4 = m
  a11, a1 = a

  return sqrt((G * (a1 * m1 * m2 + a11 * (m1 + m2) * m3) * (m1 + m2 + m3 +
    m4)) / (a1 * a11 * (m1 + m2 + m3) * m4))

def vcrit_tripbin(m, a):
  '''
  Calculate the critical velocity of a binary scattering off of a
  hierarchical triple.  (The critical velocity is the velocity at which the
  total energy of the system is zero.)

  Inputs:
    m -- a tuple containing the five masses of the stars as astropy
      quantities in the following order:
        (star 1 of the inner binary, star 2 of the inner binary, tertiary,
        star 1 of the interloping binary, star 2 of the interloping binary)

    a -- a tuple containing the three semi-major axes of the hierarchical
      triple and the interloping binary as an astropy quantity in the
      following order: (inner semi-major axis, outer semi-major axis,
      interloping binary semi-major axis)

    Output:
      The critical velocity as an astropy quantity.
  '''

  from astropy.constants import G

  # Unpack the inputs
  m000, m001, m01, m10, m11 = m
  a00, a0, a1 = a

  return sqrt((G * (m000 + m001 + m01 + m10 + m11) * (a00 * a1 * 
    (m000 + m001) * m01 + a0 * (a1 * m000 * m001 + a00 * m10 * m11))) /
    (a0 * a00 * a1 * (m000 + m001 + m01) * (m10 + m11)))

def mardling(q_out, e_out, inc=0):
  '''Calculate the ratio between the outer periapsis distance to the inner
  semi-major axis for which the triple is stable according to the Mardling
  stability criterion.  See Mardling & Aarseth (2001) for more details,
  specifically Eq. 90.

  Inputs:
    a_in -- The inner semi-major axis

    q_out -- The mass ratio between the inner and outer binaries.
      Equivalent to m3 / (m1 + m2)

    e_out -- Eccentricity of the outer orbit.

    inc [deg ] -- Inclination of the outer binary.  

  Output:
    The critical ratio between the outer periapsis distance to the inner
    semi-major axis for which the triple is Mardling stable.
  '''

  from math import sqrt
  
  # An empirically determined constant
  C = 2.8

  ratio = C * ((1 + q_out) * (1 + e_out) / sqrt(1 - e_out))**(2./5)

  # Correct for inclination (see p. 414)
  ratio *= (1 - 0.3 * inc / 180)

  return ratio

def mardling_ecc(alpha, q_out, inc=0):
  '''Calculate the critical outer eccentricity at which a hierarchical
  triple is Mardling unstable for a given semi-major axis ratio.

  Inputs:
    alpha -- The semi-major axis ratio

    q_out -- The mass ratio between the inner and outer binaries.
      Equivalent to m3 / (m1 + m2)

    inc [deg ] -- Inclination of the outer binary.  

  Output:
    The maximum eccentricity for which the triple is stable.
  '''

  import sys
  from scipy.optimize import brentq

  f = lambda ecc: mardling(q_out, ecc, inc) - (1 - ecc) * alpha
  return brentq(f, 0, 1 - sys.float_info.epsilon)
