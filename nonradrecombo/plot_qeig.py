#!/usr/bin/env python3
# -*- coding=utf-8 -*-

import sys
import os
import pickle

import argparse
import numpy as np

import matplotlib.pyplot as plt


E_MIN = -3 
E_MAX = 3 
dq = 0.01

def plot_eig(q, eigval, q0=0, e_vbm=0, marker='o'):
   
    i_kpt = 0

    for spin, eigval_s in eigval.items():
        c = '#2166ac' if float(spin) == 1 else '#b2182b'
        c = '#666666' if float(spin) == 1 else '#b2182b'
        for i_kpt in range(len(eigval_s)):
            eigval_sk = eigval[spin][i_kpt]
            occ = eigval_sk[:,1]
            y = eigval_sk[:,0]
            x = q+float(spin)*dq - q0

            # occpied
            indx_occ = np.where(occ > 0.9)
            plt.plot([x]*len(y[indx_occ]), y[indx_occ] - e_vbm, color=c, lw=0, marker=marker)
    
            # unoccpied
            indx_unocc = np.where(occ < 0.3)
            plt.plot([x]*len(y[indx_unocc]), y[indx_unocc] - e_vbm, color=c, lw=0, marker=marker, fillstyle='none')
            
            # halfoccpied
            indx_hocc = (0.3 <= occ) & (occ <= 0.9)
            if len(indx_hocc) > 0:
                plt.plot([x]*len(y[indx_hocc]), y[indx_hocc] - e_vbm, color=c, lw=0, marker=marker, fillstyle='bottom')


def plot_eigs(qs, eigvals, q0=0, e_vbm=0):   
    for q, eigval in zip(qs, eigvals):
        plot_eig(q, eigval, q0, e_vbm)


def load(filename):
    with open(filename, 'rb') as handle:
        qs, eigvals = pickle.load(handle)
    return qs, eigvals

def main(paths, q0, e_vbm, e_min, e_max):
    '''
    '''
    qs, eigvals = load("q_eig.pickle")
    plt.ylim((e_min, e_max))
    plot_eigs(qs, eigvals, q0, e_vbm)
    plt.show()


if __name__ == '__main__':
    '''
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-e","--energy", nargs=2, type=float,
                        help="plot energy range ",default=[2, 5])
    parser.add_argument("-v","--vbm", type=float,
                        help="vbm               ",default=0.)
    parser.add_argument("-q","--q0", type=float,
                        help="Q0                ",default=0.)
    #
    args = parser.parse_args()
    e_min, e_max = args.energy
    e_vbm = args.vbm
    q0    = args.q0 

    main("q_eig.pickle", q0, e_vbm, e_min, e_max)

