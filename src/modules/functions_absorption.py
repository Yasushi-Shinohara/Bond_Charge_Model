import numpy as np
from modules.constants import *

def E2Eabs(nbond, bonddict,E,npower):
    Eabs = 0.0
    for ibond in range(nbond):
        value = bonddict[str(ibond)]
        dcart = value[9]
        Edotdcart = np.abs(np.dot(E, dcart)/np.linalg.norm(dcart))
        Eabs += (Edotdcart)**npower
    return Eabs
def thetaphi2E(theta, phi):
    E = np.zeros(3, dtype='float64')
    E[0] = np.sin(theta)*np.cos(phi)
    E[1] = np.sin(theta)*np.sin(phi)
    E[2] = np.cos(theta)
    return E
def calc_Eabs_theta_phi(nbond, bonddict, npower):
    Nphi = 19
    Ntheta = 19
    theta = np.linspace(0.0,0.25,Ntheta,dtype='float64')
    phi = np.linspace(0.0,0.25,Nphi,dtype='float64')
    theta = theta*tpi
    phi = phi*tpi
    Eabstheta = np.zeros(Ntheta,dtype='float64')
    Eabsphi = np.zeros(Nphi,dtype='float64')
    for itheta in range(Ntheta):
        E = thetaphi2E(theta[itheta],phi[0])
        Eabstheta[itheta] = E2Eabs(nbond, bonddict,E,npower)
    for iphi in range(Nphi):
        E = thetaphi2E(theta[Ntheta-1],phi[iphi])
        Eabsphi[iphi] = E2Eabs(nbond,bonddict,E,npower)
    return theta, Eabstheta, phi, Eabsphi
