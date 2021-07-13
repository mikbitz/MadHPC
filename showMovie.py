#import matplotlib
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
from matplotlib import animation
import netCDF4
#import pandas as pd
import numpy as np
from scipy import ndimage
fig=plt.figure()
fig.set_size_inches(8.5, 3.5, forward=True)
plt.subplots_adjust(left=0.1, right=0.9, bottom=0., top=1.0)
plt.show()
base_dir="/home/moke/working/"
varName="totalInfected"
runNumber="119"
experiment="experiment.testCovidUK"
#fp="/home/moke/working/repastHPC/mad_modified/MadingleyData-master/0.5deg/Marine/Observed/eastwards_velocity.nc4"
#fp1="/home/moke/working/repastHPC/mad_modified/output/experiment.testingv0.22/run_013/totalCohortBiomass.nc"

#fp1="/home/mb425/repastHPC/mad_modified/output/experiment.testingv0.2withoutInteraction/run_002/"+varName+".nc"
fp1=base_dir+"covid-master/output/"+experiment+"/run_"+runNumber+"/"+varName+".nc"

#fp2="/home/moke/working/Madingley/Madv89/MadingleyCPP/output/2019-03-12_15-24-51/MonthlyGridOutputs.nc"
nc1 = netCDF4.Dataset(fp1)
#nc2 = netCDF4.Dataset(fp2)
#fig, (ax1,ax2)=plt.subplots(2,1)
t =399
n=1
k=nc1[varName][0,:,:]
v=np.mean(nc1[varName][t*n,:,:])*2;
im=np.flipud(k)

img=plt.imshow(im,vmin=0,vmax=2,cmap='viridis')
def animate(i):
    print(i)
    k=nc1[varName][i*n,:,:]
    #tph=np.ones((7,7))/49
    #tph=np.array([[0.25, 0.5, 0.25], [0.5, 1, 0.5], [0.25, 0.5, 0.25]])/4.
    im=np.flipud(k)
    #im=ndimage.convolve(np.flipud(k), tph, mode='reflect')
    img.set_data(np.roll(im,0,axis=1))
    fig.canvas.draw()
    fig.canvas.flush_events()
    return img
anim = animation.FuncAnimation(fig, animate, frames=t, interval=50)
# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Total Stock Biomass'), bitrate=1800)
#anim.save('UKStockBiomass.mp4', writer=writer)
#for f in range(1,1200):
#    print(f)
#    k=nc1['totalCohortBiomass'][f,:,:]
    #k=nc1['totalCohortBiomass'][t-n,:,:]/n
    #for i in range(t-n+1,t):
    #k=k+nc1['totalCohortBiomass'][i,:,:]/n
    # 3x3 top hat smooth - note the coasts are currently a problem here!
    #tph=np.ones((3,3))/9
    #im=ndimage.convolve(np.flipud(k), tph, mode='reflect')

    #plt.clf()
    
    #plt.contourf(im,levels=[10000,20000,40000,80000,160000,320000,640000])#np.arange(0,10,0.5))
    #avg=np.mean(k)
    #med=np.median(k)
    #print("Avg per sq km %e"% avg)
    #print("Median per sq km %e"% med)
    #print("Total biomass %e"% (avg*4*3.14*6371*6371))
    #ax2.imshow(np.flipud(nc2['AbundanceDensity'][0,:,:]),vmin=0,vmax=3.e10)
    #plt.imshow(np.flipud(nc['u'][0,:,:]))
#    img.set_data(np.flipud(k))
#    fig.canvas.draw()
#    fig.canvas.flush_events()
#    pad="000"
#    if (f>9):pad="00"
#    if (f>99):pad="0"
#    if (f>999):pad=""
    #fig.savefig("frame"+pad+str(f)+".jpg",dpi=300)
#os.system("ffmpeg -r 1 -i img%01d.png -vcodec mpeg4 -y movie.mp4")
