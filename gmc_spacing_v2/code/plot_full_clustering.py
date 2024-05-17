
### SET PARAMETERS FOR PLOTTING UNCERTAINTIES ###


med = np.zeros(len(bins))
for i in range(len(bins)):
    med[i] = np.nanmedian(clust_amp[:,i][clust_amp[:,i]!=0])
    
low = np.zeros(len(bins))
for i in range(len(bins)):
    if len(clust_amp[:,i][clust_amp[:,i]!=0]) < 2:
        low[i] = 0
    else:
        low[i] = np.nanpercentile(clust_amp[:,i][clust_amp[:,i]!=0], 16)

        
high = np.zeros(len(bins))
for i in range(len(bins)):
    if len(clust_amp[:,i][clust_amp[:,i]!=0]) < 2:
        high[i] = 0
    else:
        high[i] = np.nanpercentile(clust_amp[:,i][clust_amp[:,i]!=0], 84)
        
        
        
med_cov = np.zeros(len(bins))
for i in range(len(bins)):
    med_cov[i] = np.nanmedian(clust_amp_cov[:,i][clust_amp_cov[:,i]!=0])
    
low_cov = np.zeros(len(bins))
for i in range(len(bins)):
    if len(clust_amp_cov[:,i][clust_amp_cov[:,i]!=0]) < 2:
        low_cov[i] = 0
    else:
        low_cov[i] = np.nanpercentile(clust_amp_cov[:,i][clust_amp_cov[:,i]!=0], 16)

        
high_cov = np.zeros(len(bins))
for i in range(len(bins)):
    if len(clust_amp_cov[:,i][clust_amp_cov[:,i]!=0]) < 2:
        high_cov[i] = 0
    else:
        high_cov[i] = np.nanpercentile(clust_amp_cov[:,i][clust_amp_cov[:,i]!=0], 84)

yerr = [abs(low-med), (high-med)]
yerr_cov = [abs(low_cov-med_cov), (high_cov-med_cov)]

### MAKE CLUSTERING PLOT ###
plt.figure(figsize=(8,5), dpi=125)

plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams.update({'font.family': 'serif'})
plt.rcParams.update({'font.size': 12})

fig, ax = plt.subplots(1,2, figsize=(10,4), dpi=150)
ax[0].axvspan(158, 505, alpha=.3, hatch='//', color='gray', zorder=1)
ax[0].axvspan(158, 505, alpha=.3, hatch='//', color='None', edgecolor='gray', zorder=1)
ax[0].plot(bins, med_cov, label='Flat', c='navy', zorder=5, lw=3)
ax[0].plot(bins, med, label='Exponential', c='red', zorder=5, lw=3)
ax[0].fill_between(bins, low, high, color='red', edgecolor='None', alpha=0.2, zorder=3)

ax[0].vlines(305, ymin=-1, ymax=5, ls=':', color='black', alpha=0.7)

#ax[0].errorbar(bins, med, yerr=yerr, ls='none', capsize=3, c='red', zorder=4)
ax[0].hlines(y=1, xmin=-1000, xmax=7000, color='black', ls='--', lw=2, alpha=0.6, zorder=2)
ax[0].set_xlim(0,1300)
ax[0].tick_params(direction='in', which='both')
ax[0].set_xlabel('Cloud Spacing (pc)')
ax[0].set_ylabel(r'1 + $\omega (\Delta r)$')
ax[0].grid(axis='both', zorder=0, alpha=1, ls=':')
ax[0].set_ylim(0.25,4)
ax[0].legend()

ax[1].vlines(305, ymin=-3, ymax=15, ls=':', color='black', alpha=0.7)
ax[1].plot(bins, med_cov, label='Flat', c='navy', zorder=5, lw=3)
ax[1].plot(bins, med, label='Exponential', c='red', zorder=5, lw=3)
ax[1].fill_between(bins, low_cov, high_cov, color='navy', edgecolor='None', alpha=0.1, zorder=3)
ax[1].hlines(y=1, xmin=-1000, xmax=7000, color='black', ls='--', lw=2, alpha=0.6, zorder=2)
ax[1].set_xlim(0,6300)
ax[1].set_ylim(0,6)
ax[1].tick_params(direction='in', which='both')
ax[1].set_xlabel('Cloud Spacing (pc)')
ax[1].set_ylabel(r'1 + $\omega (\Delta r)$')
ax[1].grid(axis='both', zorder=0, alpha=1, ls=':')
ax[1].legend()

plt.savefig(home_dir + 'gmc_spacing/figures/clustering_amp_2panel.pdf', bbox_inches='tight')


fig = plt.figure(figsize=(5,4), dpi=150)
plt.plot(bins, med_cov, label='Flat', c='navy', zorder=5, lw=3)
plt.plot(bins, med, label='Exponential', c='red', zorder=5, lw=3)
plt.fill_between(bins, low, high, color='red', edgecolor='None', alpha=0.2, zorder=3)
plt.vlines(305, ymin=-1, ymax=5, ls=':', color='black', alpha=0.7)
#ax[0].errorbar(bins, med, yerr=yerr, ls='none', capsize=3, c='red', zorder=4)
plt.hlines(y=1, xmin=-1000, xmax=7000, color='black', ls='--', lw=2, alpha=0.6, zorder=2)
plt.xlim(0,6300)
plt.tick_params(direction='in', which='both')
plt.xlabel('Cloud Spacing (pc)')
plt.ylabel(r'1 + $\omega (\Delta r)$')
plt.grid(axis='both', zorder=0, alpha=1, ls=':')
plt.ylim(0.6,1.4)
plt.legend()
plt.savefig(home_dir + 'gmc_spacing/figures/clustering_stability.pdf', bbox_inches='tight')
