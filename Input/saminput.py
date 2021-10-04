#!/usr/bin/env python
# coding: utf-8
import numpy as np
import xarray as xr
def gen_omega(z, t, q, pres_sfc, D, ztop):
    #----Description of the function------
    # Input:
    # z: height (m)
    # t: potential temperature (K)
    # q: specific humidity (kg/kg)
    # pres_sfc: surface pressure (hPa)
    # D: divergence (s-1)
    # ztop: height above which the vertical velocity is a constant (m)

    #generate empty arrays for holding 
    nzm = len(z)
    presr = np.zeros(nzm)
    presi = np.zeros(nzm)
    pres = np.zeros(nzm)
    tv = np.zeros(nzm)
    zi = np.zeros(nzm)
    rho = np.zeros(nzm)

    #compute heights between original levels: zi
    dz = np.diff(np.insert(z, 0, 0))
    dz = 0.5*(dz + np.roll(dz, -1))
    dz[0] = dz[0] + 0.5*z[0]
    dz[nzm-1] = z[nzm-1]-z[nzm-2]

    for k in range(1, nzm):
        zi[k] = zi[k-1] + dz[k-1]

    #compute pressure
    rgas = 287.
    cp = 1005.
    ggr = 9.81 
    epsv = 0.61

    presr[0]=(pres_sfc/1000.)**(rgas/cp)
    presi[0]=pres_sfc

    for k in range(0, nzm-1):
        q[k]=q[k]*1e-3
        tv[k]=t[k]*(1.+epsv*q[k])
        presr[k+1]=presr[k]-ggr/cp/tv[k]*(zi[k+1]-zi[k])
        presi[k+1]=1000.*presr[k+1]**(cp/rgas)
        pres[k] = np.exp(np.log(presi[k])+np.log(presi[k+1]/presi[k])*(z[k]-zi[k])/(zi[k+1]-zi[k]))

    for k in range(0, nzm-1):
         rho[k] = (presi[k]-presi[k+1])/(zi[k+1]-zi[k])/ggr*100.

    pres[-1] = presi[-1]

    #generate vertical velocity (m/s) and omega (Pa/s) 
    w = np.zeros(nzm)
    w[z <=ztop] = -D*z[z <=ztop]
    w[z > ztop] = -D*ztop
    omega = -w*rho*ggr

    #concatenate atmosphere until pres = 0 pa
    pres_conca = np.linspace(pres[-1], 0., 1000)
    pres_full = np.concatenate((pres[:-2], pres_conca)) 

    #Get the omega of the concatenated layer, which has no observations. 
    #I just extroplolate the omega. This is not important since the top of 
    #model domain is typically lower than the top of the observation
    ind = np.where(z > ztop)[0]
    coef_lin = np.polyfit(pres[ind], omega[ind], 1)
    omega_conca = omega[-2] + coef_lin[0]*(pres_conca - pres[-2])
    omega_full = np.concatenate((omega[:-2], omega_conca)) 

    pres_full = np.concatenate(([pres_sfc], pres_full)) 
    omega_full = np.concatenate(([0], omega_full))

    #add pres coordinate to omega to make it a dataArray type
    omega_out = xr.DataArray(np.flip(omega_full),
                            dims = "lev",
                            coords = {"lev": np.flip(100.*pres_full)})

    return omega_out