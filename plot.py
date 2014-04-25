#! /usr/bin/env python

#
# Routines to help with plotting
#

def sample(X, Y, n_points=1000, xlog=False, ylog=False, bins=10):
  '''
  Sample data points to plot.  The function will sample over the full range
  of both X and Y so that no features are lost

  Inputs:
    X -- A numpy array containing the x-coordinates
    Y -- A numpy array containing the y-coordinates
    
  Optional inputs:
    npoints (optional) -- The number of data points to sample
    xlog -- Sample the x-axis linearly if False, logarithmically if True
    ylog -- Sample the y-axis linearly if False, logarithmically if True
    bins -- The number of bins to sample in

  Outputs:
    X_samp -- A numpy array containing the sampled x-coordinates
    Y_samp -- A numpy array containing the sampled y-coordinates
  '''

  import numpy as np
  import random as rand

  # Make sure that the two lists have the same length:
  if len(X) != len(Y):
    raise ValueError('plot.sample(): X and Y must have same length')

  if len(X) <= n_points:
    # If the data set is less than or equal to the number of data points
    # we're sampling, we can return the data set unchanged.
    X_samp = X
    Y_samp = Y
  else:
    # Determine if we're splitting the data linearly or logarithmically
    if xlog:
      xsplit_func = np.logspace
    else:
      xsplit_func = np.linspace

    if ylog:
      ysplit_func = np.logspace
    else:
      ysplit_func = np.linspace

    X_samp = np.zeros(n_points)
    Y_samp = np.zeros(n_points)

    xintervals = xsplit_func(min(X), max(X), num=bins, endpoint=False)
    yintervals = ysplit_func(min(Y), max(Y), num=bins, endpoint=False)

    for i, elem in enumerate(xintervals):
      binX = X[(elem < X) & (X < xintervals[i+1])]
      binY = Y[(elem < X) & (X < xintervals[i+1])]

      n_samp = n_points / 2 / bins
      rand_indices = rand.sample(range(len(binX)), n_samp)
      binX_samp = binX[rand_indices]
      binY_samp = binY[rand_indices]

      X_samp[i*n_samp:(i+1)*n_samp] = binX_samp
      Y_samp[i*n_samp:(i+1)*n_samp] = binY_samp

    for i, elem in enumerate(yintervals):
      binX = X[(elem < Y) & (Y < yintervals[i+1])]
      binY = Y[(elem < Y) & (Y < yintervals[i+1])]

      n_samp = n_points / 2 / bins
      rand_indices = rand.sample(range(len(binY)), n_samp)
      binX_samp = binX[rand_indices]
      binY_samp = binY[rand_indices]

      X_samp[nbins*n_samp + i*n_samp:nbins*n_samp + (i+1)*n_samp] = binX_samp
      Y_samp[nbins*n_samp + i*n_samp:nbins*n_samp + (i+1)*n_samp] = binY_samp

    if len(X_samp) < n_points:
      rand_indices = rand.sample(X, n_points - 2*nbins*n_samp)
      X_samp[2*nbins*n_samp+1:] = X[rand_indices]
      Y_samp[2*nbins*n_samp+1:] = Y[rand_indices]

    # Now sort X and Y along the X axis.
    sort_i = np.argsort(X_samp)
    X_samp = X_samp[sort_i]
    Y_samp = Y_samp[sort_i]
  
  return X_samp, Y_samp
