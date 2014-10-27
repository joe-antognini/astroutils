#! /usr/bin/env python

#
# stats
#
# A collection of statistics utilities
#

def wilson_score(p, n, z=1.):
  '''Calculate the Wilson score interval to estimate the uncertainty on the
  estimate of the frequency of an event with a binomial distribution.

  Inputs:
    p: number between 0 and 1
      The estimated probability

    n: non-negative number
      The number of data points.  If n is not provided as an int, the
        function will round down to the nearest int.

    z: non-negative number
      The percentile of the standard normal distribution 
          (i.e., the number of sigmas)
  
  Outputs:
    A tuple which contains the lower and upper bounds.
  '''
  from math import sqrt

  # Typechecking
  if not isinstance(p, (int, float, long)):
    raise TypeError('wilson_score(): p must be a number')
  elif not isinstance(p, (int, float, long)):
    raise TypeError('wilson_score(): n must be a number')
  elif not isinstance(z, (int, float, long)):
    raise TypeError('wilson_score(): z must be a number')
  elif p > 1 or p < 0:
    raise ValueError('wilson_score(): p must be between 0 and 1')
  elif n < 0:
    raise ValueError('wilson_score(): n must be non-negative')
  elif z < 0:
    raise ValueError('wilson_score(): z must be non-negative')
  elif not isinstance(n, int):
    n = int(n)

  discriminant = z * sqrt(p / n * (1 - p) + z**2 / (4 * n**2))

  low_bound = 1 / (1 + z**2 / n) * (p + z**2 / (2 * n) - discriminant)
  up_bound = 1 / (1 + z**2 / n) * (p + z**2 / (2 * n) + discriminant)

  return (low_bound, up_bound)
