import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
import netCDF4
import pandas as pd
import numpy as np
base_dir="/home/moke/working/"
varName="totalCohortBiomass"
runNumber="318"
experiment="experiment.testInteraction"

fp1=base_dir+"repastHPC/mad_modified/output/"+experiment+"/run_"+runNumber+"/"+varName+".nc"
nc1 = netCDF4.Dataset(fp1)
k=nc1['totalCohortBiomass'][1,:,:]/1000#/60/1000 # divide by 1000 converts to tonnes
print(np.mean(k))
k=np.roll(k,0,axis=1)
#annual average over last two years (for 12 year run)
#NB range goes from 99 to 120-1
#for i in range(60,119):
#   k=k+nc1['totalCohortBiomass'][i,:,:]/60/1000
fig=plt.figure()
fig.set_size_inches(8.5, 3.5, forward=True)
#fig.suptitle("Madingley Model Total NonVegetative Biomass")
ax=fig.add_axes()

plt.subplots_adjust(left=0.1, right=0.9, bottom=0., top=1.0)

p=plt.imshow(k,origin='lower',vmin=0,vmax=1500,cmap='bwr')#,extent=[-180,180,-65,65])


#x=np.arange(-180,180,2)
#nx=x.shape[0]
#no_labels=10
#step_x=int(nx/ (no_labels -1 ))
#x_positions = np.arange(0,nx,step_x)
#x_labels=x[::step_x]
#plt.xticks=(x_positions,x_labels)
cb=plt.colorbar(p,shrink=0.5)
cb.ax.set_ylabel('Biomass (Tonne/sq.km.)')
plt.tight_layout()
plt.show()
#fig.savefig("TotalBiomass.jpg",dpi=300)
