# -*- coding: cp1252 -*-
################################   QuakespapeDB   ##########################

from scipy import *
from numpy import *
import linecache
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import os, sys
import re
from scipy.signal import lfilter
#from signal_utilities import BUTTERWORTH
from numpy.fft import fft, fftfreq, fftshift,rfft,irfft,rfftfreq
#from tompy import time_history_plot
import scipy.ndimage
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import interpolate
import random
import bisect
import matplotlib.patches as patches
import math
import collections
from matplotlib.backend_bases import MouseEvent
################################ work directory ##############################
from os import walk
import glob
SIGN= glob.glob("DHrec/*.csv")

################ Instruction #####################
L=2.50 #Distace between shot and bore-hole
Sml=4#rigidity of the smooth coherent line in the intercept velocity (max=4, min=0)
################################################################################
###################### Import picking #######################################
#-----import time Z
Pz=genfromtxt("OUTPUT/pick_Z.asc",usecols=0)
Pz=insert(Pz, 0, 0.)
Tz=genfromtxt("OUTPUT/pick_Z.asc",usecols=1)
Tz=insert(Tz, 0, 0.)

#-----import time Y
Py=genfromtxt("OUTPUT/pick_Y.asc",usecols=0)
Py=insert(Py, 0, 0.)
Ty=genfromtxt("OUTPUT/pick_Y.asc",usecols=1)
Ty=insert(Ty, 0, 0.)

#-----import time X
Px=genfromtxt("OUTPUT/pick_X.asc",usecols=0)
Px=insert(Px, 0, 0.)
Tx=genfromtxt("OUTPUT/pick_X.asc",usecols=1)
Tx=insert(Tx, 0, 0.)

#************ P-picking correction ******************************
Hh=[Pz,Tz]
Hh=zip(*Hh)
Hh.sort(key=lambda x: x[0])
Hh=zip(*Hh)
Hh=array(Hh, float)

#******************************
alpha= arctan(Pz/L)
Thc=Tz*sin(alpha)
Phc=Pz
Hhc=[Phc,Thc]

Hhc=zip(*Hhc)
Hhc.sort(key=lambda x: x[0])
Hhc=zip(*Hhc)
Hhc=array(Hhc, float)

######mediane###############

counter=collections.Counter(Hhc[0])
#print(counter)

AA=(counter.values())

PzMp=[]
TzMp=[]
for i in range (1, len(AA)+1):
    print sum(AA[:i-1])
    P=Hhc[0][sum(AA[:i-1]):sum(AA[:i])]
    PzMp.append(median(P))
    T=Hhc[1][sum(AA[:i-1]):sum(AA[:i])]
    TzMp.append(median(T)) 

#************ P-picking:1D Coherent Interpolation ***********************

NN=len(TzMp)

f0 = InterpolatedUnivariateSpline(PzMp, TzMp, k=1)      
xnew = linspace(0, NN, num=2*NN, endpoint=True)

f1 = InterpolatedUnivariateSpline(xnew, f0(xnew), k=1)      
xnew = linspace(0, NN, num=int(NN/Sml), endpoint=True)

f2 = InterpolatedUnivariateSpline(xnew, f1(xnew), k=1)      
xnew = linspace(0, NN, num=NN+1, endpoint=True)


Vint_p=xnew[1]/ (f2(xnew)[1:]-f2(xnew)[:-1])
print Vint_p

#************ P-picking:PLOT ***********************
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1.plot(Tz, Pz,'bs', label='Picking')
ax1.set_xlabel('P-Wave Travel Time (s)')
ax1.set_ylabel('Depth (m)')
ax1.grid(True)
ax1.set_xlim(xmin=0.)
ax2.plot(Thc, Phc,'rs', label='P Corrected traveltime')
ax2.plot(TzMp, PzMp, "yo",label=' Median')
ax2.plot(f2(xnew), xnew, "--", linewidth=1.5, label='Smoothed coherent line')
ax2.set_xlabel('P-Wave Travel Time (s)')
ax2.set_ylabel('Depth (m)')
ax2.grid(True)
ax2.set_xlim(xmin=0.)
legend = ax1.legend(loc='lower left', shadow=True, fontsize=13)
legend = ax2.legend(loc='lower left', shadow=True, fontsize=13)
plt.gca().invert_yaxis()
f.canvas.set_window_title("P-wave Traveltime Analysis")

plt.show()
f.savefig('OUTPUT\P-Waves picking-PLOT.png')
Smot_Tp=f2(xnew)
Smot_Pp=xnew

########################################################################################################
#                              S-WAVES
########################################################################################################

#************ S-picking correction ******************************

Ph=concatenate((Px,Py))
Th=concatenate((Tx,Ty))
Hh=[Ph,Th]

Hh=zip(*Hh)
Hh.sort(key=lambda x: x[0])
Hh=zip(*Hh)
Hh=array(Hh, float)

#******************************
alpha= arctan(Ph/L)
Thc=Th*sin(alpha)
Thcc=Th*1.
Phc=Ph
Hhc=[Phc,Thc]
Hhcc=[Phc,Thcc]
#***************************
Hhc=zip(*Hhc)
Hhc.sort(key=lambda x: x[0])
Hhc=zip(*Hhc)
Hhc=array(Hhc, float)

#*************riporto primi arrivi S per damping (senza correzione tempi)**************
Hhcc=zip(*Hhcc)
Hhcc.sort(key=lambda x: x[0])
Hhcc=zip(*Hhcc)
Hhcc=array(Hhcc, float)
#******************************************************************************************

######mediane###############

counter=collections.Counter(Hhc[0])
#print(counter)

AA=(counter.values())

PhM=[]
ThM=[]
ThMuc=[]
for i in range (1, len(AA)+1):
    print sum(AA[:i-1])
    P=Hhc[0][sum(AA[:i-1]):sum(AA[:i])]
    PhM.append(median(P))
    T=Hhc[1][sum(AA[:i-1]):sum(AA[:i])]
    ThM.append(median(T))
    #print (P,T)

#*************riporto primi arrivi S per damping (senza correzione tempi)**************
    Tuc=Hhcc[1][sum(AA[:i-1]):sum(AA[:i])]
    ThMuc.append(median(Tuc))
#print transpose((PhM[1::],ThMuc[1::]))
savetxt("OUTPUT\Vs_median_picking.asc", transpose((PhM[1::],ThMuc[1::])), delimiter=",", comments="",
        header=' Depth up to(m), median(s)', fmt='%1.5f') 
#******************************************************************************************



#************ S-picking:1D Coherent Interpolation *********************** 

NN=len(ThM)

f0 = InterpolatedUnivariateSpline(PhM, ThM, k=1)      
xnew = linspace(0, NN, num=2*NN, endpoint=True)

f1 = InterpolatedUnivariateSpline(xnew, f0(xnew), k=1)      
xnew = linspace(0, NN, num=int(NN/Sml), endpoint=True)

f2 = InterpolatedUnivariateSpline(xnew, f1(xnew), k=1)      
xnew = linspace(0, NN, num=NN+1, endpoint=True)


Vint_s=xnew[1]/ (f2(xnew)[1:]-f2(xnew)[:-1])
print Vint_s


#************ S-picking:PLOT ***********************
fs, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1.plot(Ty, Py,'rs', label='Picking on original S-waves')
ax1.plot(Tx, Px,'bs', label='Picking on treated S-waves')
ax1.set_xlabel('S-Wave Travel Time (s)')
ax1.set_ylabel('Depth (m)')
ax1.grid(True)
ax1.set_xlim(xmin=0.)
ax2.plot(Thc,Phc,'rs',label='S Corrected traveltime')
ax2.plot(ThM, PhM, "yo",label='Median')
ax2.plot(f2(xnew), xnew, "--", linewidth=1.5, label='Smoothed coherent line')
ax2.set_xlabel('S-Wave Travel Time (s)')
ax2.set_ylabel('Depth (m)')
ax2.grid(True)
ax2.set_xlim(xmin=0.)
legend = ax1.legend(loc='lower left', shadow=True, fontsize=13)
legend = ax2.legend(loc='lower left', shadow=True, fontsize=13)
plt.gca().invert_yaxis()
fs.canvas.set_window_title("S-wave Traveltime Analysis")

plt.show()
fs.savefig('OUTPUT\S-Waves picking-PLOT.png')

#************ Vp and Vs interval ***********************

list1s = Vint_s
list2s = Vint_s
results = [None]*(len(list1s)+len(list2s))
results[::2] = list1s
results[1::2] = list2s

list1p = Vint_p
list2p = Vint_p
resultp = [None]*(len(list1p)+len(list2p))
resultp[::2] = list1p
resultp[1::2] = list2p



list1t = xnew[:-1]
list2t = xnew[1:] 
result = [None]*(len(list1t)+len(list2t))
result[::2] = list1t
result[1::2] = list2t




fvp, ax1=plt.subplots()
ax1.plot(resultp,result, linewidth=2.5, label='Vp_log by Smoothed coherent line')
ax1.set_xlabel('P-Wave velocity (m/s)')
ax1.set_ylabel('Depth (m)')

ax1.plot(results,result, linewidth=2.5, label='Vs-log by Smoothed coherent line', color='red')
ax1.set_xlabel('S-Wave velocity (m/s)')
ax1.set_ylabel('Depth (m)')

major_ticks_y_vs = arange(0, 4500, 500)
minor_ticks_y_vs = arange(0, 4500, 100) 

major_ticks_x_z = arange(0, len(xnew[1:]+1), 1)
minor_ticks_x_z = arange(0, len(xnew[1:]+1), 1)                                               

ax1.set_xticks(major_ticks_y_vs)                                                       
ax1.set_xticks(minor_ticks_y_vs, minor=True)                                           
ax1.set_yticks(major_ticks_x_z )
ax1.tick_params(axis='both', which='major', labelsize=10)
#ax1.set_yticks(minor_ticks_x_z, minor=True)
#ax1.set_xticks(fontsize=12, rotation=90)

ax1.grid(True)
#ax1.set_xscale("log", nonposx='clip')
#ax1.set_xlim(0, 2000)
#ax1.xaxis.set_ticks=2

                                         

legend = ax1.legend(loc='upper right', shadow=True, fontsize=11)
plt.gca().invert_yaxis()
plt.show()

fvp.savefig('OUTPUT\Interval VP-VS.png')
Smot_Ts=f2(xnew)
Smot_Ps=xnew


Poiss= (Vint_p**2-(2*Vint_s**2))/(2*(Vint_p**2-Vint_s**2))

savetxt("OUTPUT\Vs_VP_coherent_interval.asc", transpose((xnew[1:-1],Vint_p[:-1],Vint_s[:-1],Poiss[:-1])), delimiter=",         ", comments="",
        header=' Depth up to(m),  VP-wave (m/s), VS-wave (m/s), Poisson (--)', fmt='%1.2f') 


################ Picking domocrome VS ###########################################

TTs=[]
PPs=[]
TTs.append(0.)
PPs.append(0.)
pickS=[]

fig2, ax2 = plt.subplots()
ax2.plot(ThM, PhM, "yo",label='Median of picking')
ax2.plot(Smot_Ts, Smot_Ps, "--", linewidth=1.0, label='Smoothed coherent line')
ax2.set_xlabel('S-Wave Travel Time (s)')
ax2.set_ylabel('Depth (m)')
ax2.set_xlim(xmin=0.)
#plt.gca().invert_yaxis()
legend = ax2.legend(loc='lower left', shadow=True, fontsize=13)
plt.grid()

plt.draw()

line1, = ax2.plot(TTs,PPs,marker='+', markersize=20,color="black" )



def onclick(event):
    global iz, ix, iy
    iz,ix,iy = event.inaxes.get_title, event.xdata, event.ydata
    button=event.button
    
        
    x=event.xdata
    y=event.ydata
    new_click_x=min(TTs, key=lambda xx:abs(xx-event.xdata))
    
    
    if (abs(new_click_x-event.xdata))<0.002:
        Ax=TTs.index(new_click_x )
        
        TTs[Ax]=x       
        PPs[Ax]=y
        
        ax2=plt.gca()
        
    else:
        TTs.append(x)
        PPs.append(y)
        
    ax2=plt.gca()
    
    
    ax2.cla()
    ax2.plot(ThM, PhM, "yo",label='S Median')
    ax2.plot(Smot_Ts, Smot_Ps, "--", linewidth=1.0, label='Smoothed coherent line')
    ax2.set_xlabel('S-Wave Travel Time (s)')
    ax2.set_ylabel('Depth (m)')
    ax2.set_xlim(xmin=0.)
    
    
    line1,= ax2.plot(TTs,PPs,marker='+',markersize=20, color="black",label=' S_wave:Linear Time-to-depth envelope' )
    legend = ax2.legend(loc='lower left', shadow=True, fontsize=13)
    plt.grid()
    fig2.canvas.draw() 
    
 
connection_id = fig2.canvas.mpl_connect('button_press_event', onclick)

plt.gca().invert_yaxis()

plt.show()
fig2.savefig('OUTPUT\VS-domocrone.png')


VS_ly=[]
for i in range (0,len(PPs)-1):

    Vs_ly=(PPs[i+1]-PPs[i])/ (TTs[i+1]-TTs[i])
    VS_ly.append(Vs_ly)


################ Picking domocrome VP ###########################################
print "ok"
print pickS
print "ok"

TTs=[]
PPps=[]
TTs.append(0.)
PPps.append(0.)
pickP=[]

fig2, ax2 = plt.subplots()
ax2.plot(TzMp, PzMp, "yo",label='Median of picking')
ax2.plot(Smot_Tp, Smot_Pp, "--", linewidth=1.0, label='Smoothed coherent line')
ax2.set_xlabel('P-Wave Travel Time (s)')
ax2.set_ylabel('Depth (m)')
ax2.set_xlim(xmin=0.)
#plt.gca().invert_yaxis()
legend = ax2.legend(loc='lower left', shadow=True, fontsize=13)
plt.grid()

plt.draw()

line1, = ax2.plot(TTs,PPps,marker='+', markersize=20,color="black" )



def onclick(event):
    
    global iz, ix, iy
    iz,ix,iy = event.inaxes.get_title, event.xdata, event.ydata
    
    button=+event.button
        
    x=event.xdata
    
    print button
    new_click_x=min(TTs, key=lambda xx:abs(xx-event.xdata))
    
    
    if (abs(new_click_x-event.xdata))<0.002 :
        Ax=TTs.index(new_click_x )
        
        TTs[Ax]=x       
        #PPps[Ax]=y
        
        ax2=plt.gca()
        
    else:
        y=PPs[len(TTs)]
        TTs.append(x)
        PPps.append(y)
        
    ax2=plt.gca()
    
    
    ax2.cla()
    ax2.plot(TzMp, PzMp, "yo",label='S Median')
    ax2.plot(Smot_Tp, Smot_Pp, "--", linewidth=1.0, label='Smoothed coherent line')
    ax2.set_xlabel('P-Wave Travel Time (s)')
    ax2.set_ylabel('Depth (m)')
    ax2.set_xlim(xmin=0.)
    
    line1,= ax2.plot(TTs,PPps,marker='+',markersize=20, color="black",label=' P_wave:Linear Time-to-depth envelope' )
    legend = ax2.legend(loc='lower left', shadow=True, fontsize=13)
    plt.grid()
    fig2.canvas.draw() 
           
     
connection_id = fig2.canvas.mpl_connect('button_press_event', onclick)

plt.gca().invert_yaxis()


#plt.pause(0.001)
plt.show()
fig2.savefig('OUTPUT\VP-domocrone.png')


VP_ly=[]
for i in range (0,len(PPps)-1):

    Vp_ly=(PPps[i+1]-PPps[i])/ (TTs[i+1]-TTs[i])
    VP_ly.append(Vp_ly)

################################ Vp, Vs - Depth PLOT  #################################### 
DEPTH=[]
VP_LY=[]
VS_LY=[]
for i in range (0, len(VP_ly)):
    VP_LY.append(VP_ly[i])
    VS_LY.append(VS_ly[i])
    DEPTH.append(PPps[i])
    VP_LY.append(VP_ly[i])
    VS_LY.append(VS_ly[i])
    DEPTH.append(PPps[i+1])



fig3, ax3 = plt.subplots()
ax3.plot(VS_LY, DEPTH, linewidth=2.5, color="red", label='S-waves (Vs) ')
ax3.plot(VP_LY, DEPTH, linewidth=2.5, color="blue", label='P-waves (Vp) ')
ax3.set_xlabel('Velocity (m/s)')
ax3.set_ylabel('Depth (m)')
ax3.set_xlim(xmin=0.)
plt.gca().invert_yaxis()
legend = ax3.legend(loc='lower left', shadow=True, fontsize=13)

major_ticks_y_vs = arange(0, 4500, 500)
minor_ticks_y_vs = arange(0, 4500, 100) 

major_ticks_x_z = arange(0, len(xnew[1:]+1), 1)
minor_ticks_x_z = arange(0, len(xnew[1:]+1), 1)

ax3.set_xticks(major_ticks_y_vs)                                                       
ax3.set_xticks(minor_ticks_y_vs, minor=True)                                           
ax3.set_yticks(major_ticks_x_z )
ax3.tick_params(axis='both', which='major', labelsize=10)

plt.grid()

plt.show()

fig3.savefig('OUTPUT\VP - VS Profiles.png')
################################ Print Vp, Vs Profile  ##########################
H=[]# Thickness
for i in range (0, len(VP_ly)):
	H.append(PPps[i+1]-PPps[i])


Z=PPps[1:]#Depth

VP_ly=array(VP_ly,float)
VS_ly=array(VS_ly,float)
#POISS= ((0.5*(VP_ly/VS_ly)**2)-(VP_ly/VP_ly))/(((VP_ly/VS_ly)**2)-(VP_ly/VP_ly))#Poisson
POISS= (VP_ly**2-(2*VS_ly**2))/(2*(VP_ly**2-VS_ly**2)) #               ((0.5*(VP_ly/VS_ly)**2)-(VP_ly/VP_ly))/(((VP_ly/VS_ly)**2)-(VP_ly/VP_ly))

savetxt("OUTPUT\Vs_VP_profile.asc", transpose([Z[:-1],H[:-1],VP_ly[:-1],VS_ly[:-1],POISS[:-1]]), delimiter=",         ", comments="",
        header=' Depth up to(m), Tickness (m), VP-wave (m/s), VS-wave (m/s), Poisson (--)',
        footer='>%1.2f,     >%1.2f,       %1.2f,         %1.2f,         %1.2f'%((Z[-1]),(H[-1]),(VP_ly[-1]),(VS_ly[-1]),(POISS[-1])), fmt='%1.2f') 
















