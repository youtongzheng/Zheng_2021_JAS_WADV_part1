#!/usr/bin/env python
# coding: utf-8

#Usage: python Postprc_SAM.py <dir> <case_name>

import numpy as np
import xarray as xr
import sys
import glob
import os
from yzlcl import *

# In[2]:
# Parse inputs
if len(sys.argv) !=3:
    print('ERROR---Usage: python Postprc_SAM.py ' +
        '<dir> <case_name>' )
    exit()

stat_dir = sys.argv[1] + sys.argv[2] 
dest = stat_dir 

# In[3]:
Rd     = 287.              # Gas constant for dry air, J/kg/K
Rv     = 461.              # Gas constant for water vapor, J/kg/K
cp     = 1004.             # Specific heat of air, J/kg/K
grav    = 9.81              # Gravity acceleration, m/s2
Lv     = 2.5104e6          # Latent heat of condensation, J/kg

# Create dictionaries for holding postprocessed variables
data_vars = {}
coords = {}
attrs = {}

# Open statistics files
f = xr.open_mfdataset(stat_dir + '/*N100.nc', combine = 'by_coords', 
                         decode_cf = False)

nt = f['time'].size
nz = f['z'].size

z = f['z']
p = f['p']

# Save required coordinates
coords['z'] = f['z']
coords['time'] = f['time']

z2D = f['z'].expand_dims(time = f['time'].size)

#Determine cloud base/top
cld = f['CLD'] # cloud fraction

data_vars['zcb'] = z2D.where(cld >= 0.5).min(dim = 'z', skipna = 'True')

data_vars['zcb'].attrs = {
    'long name': 'cloud base height of Sc deck',
    'units': 'm',
    'source': 'OUT_STAT',
    'criteria': 'lowest grid with cloud fraction > 0.5'
}

data_vars['zct'] = z2D.where(cld >= 0.5).max(dim = 'z', skipna = 'True')

data_vars['zct'].attrs = {
    'long name': 'cloud top height of Sc deck',
    'units': 'm',
    'source': 'OUT_STAT',
    'criteria': 'highest grid with cloud fraction > 0.5'
}

#Determine caping inversion base/top
thl2 = f['THL2']

max_thl2 = thl2.max(dim = 'z', skipna = 'True')

data_vars['zinvt'] = z2D.where(thl2 >= 0.05*max_thl2).max(dim = 'z', skipna = 'True')
data_vars['zinvt'].attrs = {
    'long name': 'top of caping inversion',
    'units': 'm',
    'source': 'OUT_STAT',
    'criteria': 'Based on van der Dussen et al., 2014'
}

data_vars['zinvb'] = z2D.where(thl2 >= 0.05*max_thl2).min(dim = 'z', skipna = 'True')
data_vars['zinvb'].attrs = {
    'long name': 'base of caping inversion',
    'units': 'm',
    'source': 'OUT_STAT',
    'criteria': 'Based on van der Dussen et al., 2014'
}

zcb = data_vars['zcb']
zct = data_vars['zct']
zinvb = data_vars['zinvb']
zinvt = data_vars['zinvt']

#determine quantities at inversion base/top
thl  = f['THETAL']
thv  = f['THETAV']
qt   = f['QT']/1000.

trho_tmp = f['TABSOBS']*(1. + 0.61*f['QV']/1000. - f['QC']/1000.)
buo = (grav/trho_tmp)*(trho_tmp+ grav*z2D/cp)

data_vars['thl_zinvt'] = thl.sel(z = zinvt, method="nearest")
data_vars['thl_zinvb'] = thl.sel(z = zinvb, method="nearest")

data_vars['thv_zinvt'] = thv.sel(z = zinvt, method="nearest")
data_vars['thv_zinvb'] = thv.sel(z = zinvb, method="nearest")

data_vars['qt_zinvt'] = qt.sel(z = zinvt, method="nearest")
data_vars['qt_zinvb'] = qt.sel(z = zinvb, method="nearest")

data_vars['b_zinvt'] = buo.sel(z = zinvt, method="nearest")
data_vars['b_zinvb'] = buo.sel(z = zinvb, method="nearest")

# determine PBL mean
data_vars['thl_scale'] = thl.where(z2D <= data_vars['zinvb']).mean(dim='z')
data_vars['qt_scale'] = qt.where(z2D <= data_vars['zinvb']).mean(dim='z')

# determine LWP budegt
qtf  = f['QTFLUX']/Lv/f['RHO']  # Nonprecipitating water flux (Total), W/m2 -> kg/kg m/s
thlf = f['THLFLUX'] # THETAL flux (Resolved), K m/s
precf = f['PRECIP']/24/3600 # sfc precipitation m/day --> m/s
rad  = f['RADLWUP'] - f['RADLWDN'] + f['RADSWUP'] - f['RADSWDN']  # W/m2

h = data_vars['zinvt'] - data_vars['zcb']

rho_avg = f['RHO'].where(z2D >= data_vars['zcb']).where(z2D <= data_vars['zinvt']).mean(dim='z')
qv_avg = f['QV'].where(z2D >= data_vars['zcb']).where(z2D <= data_vars['zinvt']).mean(dim='z')/1000.
t_avg = f['TABSOBS'].where(z2D >= data_vars['zcb']).where(z2D <= data_vars['zinvt']).mean(dim='z')

gamma   = Lv*qv_avg/Rv/t_avg**2
eta     = (1+Lv*gamma/cp)**(-1)
Ga_ql   = grav*eta*(qv_avg/Rd/t_avg - gamma/cp)

qt_invt = qt.sel(z = zinvt, method="nearest")
qt_invb = qt.sel(z = zinvb, method="nearest")

thl_invt = thl.sel(z = zinvt, method="nearest")
thl_invb = thl.sel(z = zinvb, method="nearest")

thv_invt = thv.sel(z = zinvt, method="nearest")
thv_invb = thv.sel(z = zinvb, method="nearest")

p_invt = p.sel(z = zinvt, method="nearest")
p_invb = p.sel(z = zinvb, method="nearest")

rho_invt = f['RHO'].sel(z = zinvt, method="nearest")
rho_invb = f['RHO'].sel(z = zinvb, method="nearest")

d_qt   = qt_invt-qt_invb
d_thl  = thl_invt-thl_invb
    
exner_invt  = (p_invt/f['Ps'])**(Rd/cp)
exner_invb  = (p_invb/f['Ps'])**(Rd/cp)

qtf_ct = qtf.sel(z = zinvt, method="nearest")
qtf_cb = qtf.sel(z = zcb, method="nearest")

thlf_ct = thlf.sel(z = zinvt, method="nearest")
thlf_cb = thlf.sel(z = zcb, method="nearest")

precf_ct = precf.sel(z = zinvt, method="nearest")
precf_cb = precf.sel(z = zcb, method="nearest")

rad_ct = rad.sel(z = zinvt, method="nearest")
rad_cb = rad.sel(z = zcb, method="nearest")

### Sub ###
w_zi = f['WOBS'].sel(z = 1000.*f['ZINV'], method="nearest")
data_vars['Subs_LWP'] = (-1.0*rho_avg*h*Ga_ql*w_zi)*1000*3600      # kg/m2/s --> g/m2/h

### Ent ###
data_vars['Ent_LWP'] = rho_avg*f['we']*(eta*d_qt - gamma*eta*(exner_invt*thl_invt-exner_invb*thl_invb) - h*Ga_ql)*1000*3600

### Base ###
data_vars['Base_LWP'] = rho_avg*eta*(qtf_cb - exner_invb*gamma*thlf_cb)*1000*3600  

### Rad ###
data_vars['rad_LWP'] = eta*gamma/cp*(rad_ct - rad_cb)*1000*3600 

### Prec ###
data_vars['Prec_LWP']  = 1.0*rho_avg*(precf_ct - precf_cb)*1000*3600  

#determine other quantities
data_vars['CldRCool'] = rad_ct - rad_cb

data_vars['CldtopRCool'] = rad.where((z2D >= zcb) & (z2D <= zinvt)).max(
    dim = 'z', skipna = 'True') - rad.where((z2D >= zcb) & (z2D <= zinvt)).min(
    dim = 'z', skipna = 'True')

data_vars['kappa'] = 1. + (cp/Lv)*(d_thl/d_qt) 

data_vars['dissip_invb20'] = f['DISSIP'].where((z2D >= 0.8*zinvb) & (z2D < zinvb)).mean(dim = 'z')

data_vars['TVFLX_we'] = f['we']*(thv_invt - thv_invb)
data_vars['TVFLX_we_Wm2'] = data_vars['TVFLX_we']*(rho_invb*cp)
data_vars['Prec_wm2']  = 0.91*Lv*(precf_cb - precf_ct)

data_vars['delta_thetal']  = thl.where((z2D > 0.9*zinvb) & (z2D < zinvb)).mean(
    dim = 'z')-thl.where((z2D > 0) & (z2D < 0.1*zinvb)).mean(dim = 'z')

data_vars['delta_qt']  = (qt.where((z2D > 0.9*zinvb) & (z2D < zinvb)).mean(
    dim = 'z')-qt.where((z2D > 0) & (z2D < 0.1*zinvb)).mean(dim = 'z'))*1000.

#determine lifting condensation level (lcl) using Romps's 2014 formula
data_vars['LCL'] = yzlcl(100.*f['p'].isel(z = slice(0, 10)).mean(dim = 'z'),
    f['TABSOBS'].isel(z = slice(0, 10)).mean(dim = 'z'),
    f['QV'].isel(z = slice(0, 10)).mean(dim = 'z')/f['QSAT'].isel(z = slice(0, 10)).mean(dim = 'z')
    )
# Create output NetCDF
print('-----------------Writing postprocessed file.--------------')
xr.Dataset(data_vars, coords = coords).to_netcdf(
    dest + '/postproc_out_N100.nc')

