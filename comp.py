#!/usr/bin/python
# For license details see LICENSE.

import sys
import math
import numpy as np
import scipy.spatial.distance
import matplotlib.pyplot as plt
from scipy.stats import linregress

box_l = 2600.0

def analyzeAgglomerate(A, sigma, plot_regression=None):
    """computes fractal dimension
    assumes that max distance between any particles in agglomerate
    is less than half the box length box_l AND does not cross the boundary"""
    numPart = A.shape[0]
    # Hack to unfold the particle positions -- unnecessary; agglo.py unfolds
    # them now correctly, even if agglomerate is longer than box/2.
    # Currently, unfolding in agglo.py is deactivated.
    A = (A - (A[0,:] - box_l/2.)) % box_l

    #print("Pos after normalization:", A)

    # compute the maximal distance
    longdist = np.max(scipy.spatial.distance.pdist(A))
    # compute centre of mass
    com = np.sum(A, axis=0) / A.shape[0]
    # compute distance vector
    dist = np.apply_along_axis(np.linalg.norm, 1, A - com)

    #print("Center of mass:", com)

    diameters = []
    pcounts = []
    radogs = []
    k = 0
    rad = 0
    while k < numPart:
        rad += sigma
        distvec = dist[dist < rad]
        k = len(distvec)
        if k > 0:
            pcounts.append(k)
            diameters.append(2 * rad)
            radogs.append(np.sqrt(np.dot(distvec, distvec) / k))
    diameters = np.array(diameters)
    pcounts = np.array(pcounts)
    radogs = np.array(radogs)
    # compute Df (based on both diameter and on radius of gyration)
    Df_diam, intercept, _, _, _ = linregress(np.log(diameters), np.log(pcounts))

    if plot_regression == numPart:
        plt.figure()
        xs = np.log(diameters)
        plt.plot(xs, np.log(pcounts), "bx",
                 xs, intercept + Df_diam * xs, "g",
                 xs, intercept + 2.0 * xs, "r")
        plt.title("Agglo {}".format(numPart))

    #print("Pcounts:", pcounts)
    #print("Diams:", diameters)
    #print("LogPcounts:", np.log(pcounts))
    #print("LogRadogs:", np.log(radogs))

    radogs = radogs / sigma # Normalization by radius adapted from [Xiong and Friedlander 2001]
    Df_radog, intercept, _, _, _ = linregress(np.log(radogs), np.log(pcounts))
    return longdist, Df_diam, radogs[-1], Df_radog


A = np.loadtxt(sys.argv[1])
x = analyzeAgglomerate(A, 1.0)
print(x[3])
