#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 17:00:34 2019

@author: moke
"""

from matplotlib import pyplot as plt
import netCDF4


fig=plt.figure()
fig.set_size_inches(8.5, 3.5, forward=True)
#plt.subplots_adjust(left=0.1, right=0.9, bottom=0., top=1.0)
plt.show()

base_dir="/home/moke/working/"
varName="totalCohortBreakdown"
runNumber="309"
experiment="experiment.testInteraction"

fp1=base_dir+"repastHPC/mad_modified/output/"+experiment+"/run_"+runNumber+"/"+varName+".nc"

nc1 = netCDF4.Dataset(fp1)

for j in range(0,18): #range goes to one less than upper lim
    k=nc1['numberOfCohortsInFunctionalGroup'][0:,j]
    plt.plot(k)
    

