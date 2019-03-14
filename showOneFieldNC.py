import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
import netCDF4
import pandas as pd
import numpy as np

fp1="/home/moke/working/repastHPC/mad_modified/output/experiment.tenDegreeCheck/run_000/totalCohortBiomass.nc"
nc1 = netCDF4.Dataset(fp1)
k=nc1['totalCohortBiomass'][0,:,:]/120
print(np.mean(k))
#annual average over last two years (for 12 year run)
#NB range goes from 99 to 120-1
for i in range(1,120):
   k=k+nc1['totalCohortBiomass'][i,:,:]/120
   
p=plt.imshow(np.flipud(k),vmin=0,vmax=500000,cmap='bwr')
plt.colorbar(p)
plt.tight_layout()
plt.show()
