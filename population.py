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
  m_l = 0.01
  m_u = 150

  G_m_l = (1 + (m_l / mu)**(1 - alpha))**(1 - beta)
  G_m_u = (1 + (m_u / mu)**(1 - alpha))**(1 - beta)

  return (mu * ((u * (G_m_u - G_m_l) + G_m_l)**(1 / (1 - beta)) 
    - 1)**(1 / (1 - alpha)))
