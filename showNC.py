#import matplotlib
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
import math
import netCDF4
#import pandas as pd
import numpy as np
from scipy import ndimage
fig=plt.figure()
fig.set_size_inches(8.5, 3.5, forward=True)
plt.subplots_adjust(left=0.1, right=0.9, bottom=0., top=1.0)
plt.show()
#fp="/home/moke/working/repastHPC/mad_modified/MadingleyData-master/0.5deg/Marine/Observed/eastwards_velocity.nc4"

#fp1="/home/moke/working/repastHPC/mad_modified/output/experiment.runTests/run_252/totalCohortBiomass.nc"
fp1="/home/moke/working/repastHPC/mad_modified/output/experiment.testInteraction/run_308/totalCohortBiomass.nc"

#fp1="/home/moke/working/repastHPC/mad_modified/output/experiment.testingv0.2withInteraction/run_031/totalCohortBiomass.nc"
#fp2="/home/moke/working/Madingley/Madv89/MadingleyCPP/output/2019-03-12_15-24-51/MonthlyGridOutputs.nc"
nc1 = netCDF4.Dataset(fp1)
#nc2 = netCDF4.Dataset(fp2)
#fig, (ax1,ax2)=plt.subplots(2,1)
t = np.shape(nc1['totalCohortBiomass'])[0]-1
n=1
k=nc1['totalCohortBiomass'][t*n,:,:]/1


   
#for i in range(t-n+1,t,12):
#    k=k+nc1['totalCohortBiomass'][i,:,:]/12
# 3x3 top hat smooth - note the coasts are currently a problem here!
#tph=np.ones((3,3))/9
tph=np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])/4.
im=ndimage.convolve(np.flipud(k), tph, mode='reflect')/np.flipud(k)-1
#im=np.roll(np.flipud(k),361,axis=1)

img=plt.imshow(im,vmin=0,vmax=np.mean(im),cmap='viridis')
    #plt.clf()
    
#plt.contourf(im,levels=[10000,20000,40000,80000,160000,320000,640000])#np.arange(0,10,0.5))
avg=np.mean(k)
med=np.median(k)
rms=math.sqrt(np.mean(np.square(im)))
print("Avg per sq km %e"% avg)
print("Median per sq km %e"% med)
print("Total biomass %e"% (avg*4*3.14*6371*6371))
print("RMS error %e"% rms)
    #ax2.imshow(np.flipud(nc2['AbundanceDensity'][0,:,:]),vmin=0,vmax=3.e10)
    #plt.imshow(np.flipud(nc['u'][0,:,:]))
#    img.set_data(np.flipud(k))
#    fig.canvas.draw()
#    fig.canvas.flush_events()

