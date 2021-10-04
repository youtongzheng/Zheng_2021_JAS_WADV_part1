"""
Parameters (data paths, etc) used to make figures
"""
import matplotlib.pyplot as plt
import numpy as np

output_dir = '../SAM_RUNS'
wtg_2K = 'CWTG_SST302_NC%d'
wtg_1K = 'CWTG_SST301_NC%d'
wtg_0d5K = 'CWTG_SST300d5_NC%d'
wtg_0d25K = 'CWTG_SST300d25_NC%d'
wtg_0K = 'CWTG_SST300_NC%d'
wtg_m0d25K = 'CWTG_SST299d75_NC%d'
wtg_m0d5K = 'CWTG_SST299d5_NC%d'
wtg_m1K = 'CWTG_SST299_NC%d'
rce = 'RCE_SST300_NC%d'
rce_small = 'RCE_SST300_NC%d_SMALL'
rce_lf0 = 'RCE_SST300_NC%d_LF0'
rce_lf0_small = 'RCE_SST300_NC%d_LF0_SMALL'
rce_cc = 'RCE_SST300_NC%d_1KDAY'
rce_cc_lf0 = 'RCE_SST300_NC%d_LF0_1KDAY'
rce_wm5 = 'RCE_SST300_NC%d_WM5'
rce_wm1 = 'RCE_SST300_NC%d_WM1'
rce_wm0d5 = 'RCE_SST300_NC%d_WM0D5'
rce_wm0d25 = 'RCE_SST300_NC%d_WM0D25'
rce_w1 = 'RCE_SST300_NC%d_W1'
rce_w5 = 'RCE_SST300_NC%d_W5'
rce_w = [rce_wm5, rce_wm1, rce_wm0d5, rce_wm0d25, rce, rce_w1, rce_w5]
rce_wup = [rce, rce_w1, rce_w5]
all_sims = [wtg_2K, wtg_1K, wtg_0d5K, wtg_0d25K, wtg_0K,
           wtg_m0d25K, wtg_m0d5K, wtg_m1K, rce]
wtg_up = [wtg_0K, wtg_0d25K, wtg_0d5K, wtg_1K, wtg_2K]
wtg_down = [rce, wtg_0K, wtg_m0d25K, wtg_m0d5K, wtg_m1K]
wtg_all = [wtg_m1K, wtg_m0d5K, wtg_m0d25K, wtg_0K, wtg_0d25K, 
           wtg_0d5K, wtg_1K, wtg_2K]
sim_labels = {wtg_2K: '+2K WTG', 
              wtg_1K: '+1K WTG', 
              wtg_0d5K: '+0.5K WTG', 
              wtg_0d25K: '+0.25K WTG', 
              wtg_0K: '+0K WTG',
              wtg_m0d25K: '-0.25K WTG', 
              wtg_m0d5K: '-0.5K WTG', 
              wtg_m1K: '-1K WTG', 
              rce: 'RCE'}

def Nc_colors(n):
    return plt.get_cmap('plasma')(np.linspace(0,0.85,n))