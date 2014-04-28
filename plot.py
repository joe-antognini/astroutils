#! /usr/bin/env python

#
# Routines to help with plotting
#

def sample(X, Y, n_points=1000, return_indices=False, xlog=False, 
  ylog=False, bins=10):
  '''
  Sample data points to plot.  The function will sample over the full range
  of both X and Y so that no features are lost

  Inputs:
    X -- A numpy array containing the x-coordinates
    Y -- A numpy array containing the y-coordinates
    
  Optional inputs:
    npoints -- The number of data points to sample
    return_indices -- If true, also return the indices
    xlog -- Sample the x-axis linearly if False, logarithmically if True
    ylog -- Sample the y-axis linearly if False, logarithmically if True
    bins -- The number of bins to sample in

  Outputs:
    X_samp -- A numpy array containing the sampled x-coordinates
    Y_samp -- A numpy array containing the sampled y-coordinates
  '''

  import numpy as np
  import random as rand
  import sys

  # Make sure that the two lists have the same length:
  if len(X) != len(Y):
    raise ValueError('plot.sample(): X and Y must have same length')

  if len(X) <= n_points:
    # If the data set is less than or equal to the number of data points
    # we're sampling, we can return the data set unchanged.
    X_samp = X
    Y_samp = Y
  else:
    X_samp = np.zeros(n_points)
    Y_samp = np.zeros(n_points)
    indices = np.arange(len(X))
    indices_samp = np.zeros(n_points)

    # Sample either linearly or logarithmically in x and y
    if xlog:
      xintervals = np.logspace(np.log10(min(X)), np.log10(max(X)), num=bins)
    else:
      xintervals = np.linspace(min(X), max(X), num=bins)

    if ylog:
      yintervals = np.logspace(np.log10(min(Y)), np.log10(max(Y)), num=bins)
    else:
      yintervals = np.linspace(min(Y), max(Y), num=bins)

    total_points_sampled = 0

    # Sample points from the full range in the x-coordinate
    for i, elem in enumerate(xintervals[:-1]):
      binX = X[(elem < X) & (X < xintervals[i+1])]
      binY = Y[(elem < X) & (X < xintervals[i+1])]
      bini = indices[(elem < X) & (X < xintervals[i+1])]

      n_samp = min(n_points / 2 / bins, len(binX))
      rand_indices = rand.sample(range(len(binX)), n_samp)
      binX_samp = binX[rand_indices]
      binY_samp = binY[rand_indices]
      bini_samp = bini[rand_indices]

      # Put the sampled data points from this bin into the output array
      X_samp[total_points_sampled:total_points_sampled + n_samp] = binX_samp
      Y_samp[total_points_sampled:total_points_sampled + n_samp] = binY_samp
      indices_samp[total_points_sampled:total_points_sampled + n_samp] = \
        bini_samp
      total_points_sampled += n_samp

    # Sample points from the full range in the y-coordinate
    for i, elem in enumerate(yintervals[:-1]):
      binX = X[(elem < Y) & (Y < yintervals[i+1])]
      binY = Y[(elem < Y) & (Y < yintervals[i+1])]
      bini = indices[(elem < Y) & (Y < yintervals[i+1])]

      n_samp = min(n_points / 2 / bins, len(binX))
      rand_indices = rand.sample(range(len(binY)), n_samp)
      binX_samp = binX[rand_indices]
      binY_samp = binY[rand_indices]
      bini_samp = bini[rand_indices]

      X_samp[total_points_sampled:total_points_sampled + n_samp] = binX_samp
      Y_samp[total_points_sampled:total_points_sampled + n_samp] = binY_samp
      indices_samp[total_points_sampled:total_points_sampled + n_samp] = \
          bini_samp
      total_points_sampled += n_samp

    # If we didn't get enough points, pad the rest of the output with random
    # points
    if total_points_sampled < n_points:
      rand_indices = rand.sample(range(len(X)), n_points - total_points_sampled)
      X_samp[total_points_sampled:] = X[rand_indices]
      Y_samp[total_points_sampled:] = Y[rand_indices]
      indices_samp[total_points_sampled:] = rand_indices

    # Now sort X and Y along the X axis.
    #sort_i = np.argsort(X_samp)
    #X_samp = X_samp[sort_i]
    #Y_samp = Y_samp[sort_i]
    #if return_indices:
    #  indices_samp = indices_samp[sort_i]
  
  if return_indices:
    return X_samp, Y_samp, indices_samp
  else:
    return X_samp, Y_samp
