#! /usr/bin/env python

#
# Population library
#
# Contains routines describing stellar populations: IMFs, stellar lifetimes,
# etc.
#

def rand_imf(m_l=.01, m_u=150):
  '''
  Randomly draw a stellar mass from an IMF.  The IMF used is that
  presented in Maschberger (2013) which was designed to be an analytic
  function that matches Kroupa (2001, 2002) and Chabrier (2003).
  Specifically, this function implements Eq. (4) of Table 1 of Maschberger
  (2013). 

  Parameters:
    m_l: float, optional 
      The lower bound on the mass

    m_u: float, optional
      The upper bound on the mass

  Returns:
    m: float
      A random mass

  References:
    Chabrier, G., 2003, PASP, 115, 763
    Kroupa, P., 2001, MNRAS, 322, 231
    Kroupa, P., 2002, Sci., 295, 82
    Maschberger, T., 2013, MNRAS, 429, 1725
  '''

  import random

  u = random.random()

  # Parameters
  # (See Table 1 of Maschberger 2013)
  alpha = 2.3
  beta = 1.4
  mu = 0.2 

  G_m_l = (1 + (m_l / mu)**(1 - alpha))**(1 - beta)
  G_m_u = (1 + (m_u / mu)**(1 - alpha))**(1 - beta)

  m = (mu * ((u * (G_m_u - G_m_l) + G_m_l)**(1 / (1 - beta)) 
    - 1)**(1 / (1 - alpha)))

  return m

def stellar_lifetime(m):
  '''
  Calculate the lifetime of a star in years given its mass in solar masses.
  See Eq. (1.91) of Hansen, Kawaler, & Trimble (2004)

  Parameters:
    m: float
      The mass of the star (solar masses)
  
  Returns:
    t: float
      The main sequence lifetime of the star (years)

  References:
    Hansen, C. J., Kawaler, S. D.,, & Trimble, V. 2004, Stellar Interiors
      (New York: Springer)
  '''

  exponent = -2.9

  # At masses above m_crit, the lifetime is constant at 50 Myr.
  t_crit = 5e7
  m_crit = (t_crit / 1e10)**(1. / exponent)
  if m > m_crit:
    return t_crit
  else:
    return 1e10 * m**exponent

def wd_ifmr(mi):
  '''
  Calculate the final mass of a white dwarf given the initial mass of the
  star.  The calculation is based on Eq. (5) of Zhoa et al. (2012).  This
  equation is only valid for initial masses of 1.1 < m < 4.1.

  Parameters:
    mi: float
      Initial mass of the star (solar masses)

  Returns:
    mf: float
      Final mass of the star (solar masses)
  
  References:
    Zhoa, J.K., et al., 2012, ApJ, 746, 144
  '''

  # Make sure that the range is valid.
  if m < 1.1 or m > 4.1:
    raise ValueError
  
  mf = 0.452 + 0.073 * mi

  return mf
