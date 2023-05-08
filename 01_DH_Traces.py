# -*- coding: cp1252 -*-
################################   QuakespapeDB   ##########################

from scipy import *
from numpy import *
import linecache
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.ticker as ticker
import os, sys
import re
from scipy.signal import lfilter
#from signal_utilities import BUTTERWORTH
from numpy.fft import fft, fftfreq, fftshift,rfft,irfft,rfftfreq
#from tompy import time_history_plot
import scipy.ndimage
from scipy.interpolate import interp1d
from scipy import interpolate
import random
import bisect
from scipy import stats



################################ work directory ##############################
from os import walk
import glob
SIGN= glob.glob("DHrec/*.txt")



################################ FILTER ##############################
from scipy.signal import butter, lfilter

def butter_bandpass(lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=3):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def pos(A,value):
    A=array(A)
    i=0
    a=0
    while a<=value:
        a=A[i]
        i=i+1
    return i
    

################################ Istruzioni ##############################
DT=0.064#(msec) #simple time
dz= 1. # intergeofone distance 
cut=2.0 # cut the trace windows e.g. 2 is equal to half

################################ Iterpretazione ##############################
delay_time=0.025#sec
Vmed=580#extimate Vs (m/s) mean traveling in the hole
change=[]#change X in Y...sostituisce alcune tracce dovute a probabile rotazione..[]=nessuna sostituz.
#change=[12,13,14,15,16,17,18,19,20,21,22,23,24,25]
Inv=0# inverte l'asse X con quella Y
Zoom=[]
X_Zoom=2.0

spec=0##### inverte le tracce da fondo foro
zto3=0#### la battuta verticale diventa la terza

if Inv==1:
    A=3
    B=2
else:
    A=2
    B=3

######################################
X_Zoom=1/X_Zoom
change=array((change-ones(len(change))),int)
########## crea cartella parameter zones ################

path1 = 'OUTPUT'
if os.path.exists(path1):
    pass
else:
    os.makedirs( path1, 0755 );

path1 = 'Z_signal'
if os.path.exists(path1):
    pass
else:
    os.makedirs( path1, 0755 );


path1 = 'Y_signal'
if os.path.exists(path1):
    pass
else:
    os.makedirs( path1, 0755 );

    
path1 = 'X_signal'
if os.path.exists(path1):
    pass
else:
    os.makedirs( path1, 0755 );


################################################################################

ID_Sign=[]
Files_Sign=[]


for i,j in enumerate(SIGN):
    ID_Sign.append(i)
    Files_Sign.append(j)


#################reverse tracce#####

if spec ==1:
    ### specula l'ordine delle tracce ####
    ID_Sign=ID_Sign[::-1]
    Files_Sign=Files_Sign[::-1]

if zto3==1:
    ### considera z alla 3th battuta ####
    ID_Sign_rev=[]
    Files_Sign_rev=[]
    i=0
    
    while i<=(len(Files_Sign)-1):
        A=ID_Sign[2+i]
        ID_Sign_rev.append(A)
        B=ID_Sign[1+i]
        ID_Sign_rev.append(B)
        C=ID_Sign[0+i]
        ID_Sign_rev.append(C)

        AA=Files_Sign[2+i]
        Files_Sign_rev.append(AA)
        BB=Files_Sign[1+i]
        Files_Sign_rev.append(BB)
        CC=Files_Sign[0+i]
        Files_Sign_rev.append(CC)
        i=i+3
        
    
    ID_Sign=ID_Sign_rev
    Files_Sign=Files_Sign_rev





#DT= genfromtxt("%s"%(Files_Sign[0]),skip_header=10, delimiter=",", usecols=2)#simple time


Z=[]
Yp=[]
Yn=[]
Xp=[]
Xn=[]
Ypf=[]
Ynf=[]
Xpf=[]
Xnf=[]
Ypff=[]
Ynff=[]
Ypp=[]
Ynn=[]
Xpp=[]
Xnn=[]
lowpass=20
highpass=100
highpassp=500
fs=10000
Or=1



IZ=[]
i=2
while i<=(len(Files_Sign)-1):
    z=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=1)
    z=butter_bandpass_filter(z, lowpass, highpassp, fs, order=Or)
    z=z[:int((len(z)/cut))]
    one=ones(len(z))
    z=z-(one*(mean(z)))
    z4=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=4)
    z4=z4/max(abs(z4))
    iz=pos(abs(z4),0.02)
    IZ.append(iz)
    Z.append(z/max(abs(z)))
    print Files_Sign[i]
    
    i=i+3
################# pareggio tempi###############
aa=stats.mode(IZ)
for i in range (0, len(IZ)):
    delta_t=int(aa[0]-IZ[i])
    if delta_t>0:
        Z[i]= list([0] * delta_t)+list(Z[i])
        Z[i]=Z[i][0:-delta_t]
    elif delta_t<0:
        Z[i]=Z[i][abs(delta_t)-1:-1]
        Z[i]= list(Z[i])+list([0] *abs(delta_t))
    else:
        pass
       
    
    #print'%s Triplet Imported'%(i/3)
IYP=[]   
i=1
while i<=(len(Files_Sign)-1):
    yp=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=A)
    ypp=yp
    
    yp=butter_bandpass_filter(yp, lowpass, highpass, fs, order=Or)
    yp=yp[:int(len(yp)/cut)]
    yp=yp-(one*(mean(yp)))
    Yp.append(yp/max(abs(yp)))
    
    yp4=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=4)
    yp4=yp4/max(abs(yp4))
    iyp=pos(abs(yp4),0.02)
    IYP.append(iyp)
    #Ypf.append(butter_bandpass_filter(ypp, 70, 80, fs, order=Or))

    
    xpp=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=B)
    ypp=(ypp**2 + xpp**2)**0.5 *sign(ypp + xpp)
    ypp=ypp[:int(len(ypp)/cut)]
    ypp=ypp-(one*(mean(ypp)))
    ypp=butter_bandpass_filter(ypp, lowpass, highpass, fs, order=Or)
    Ypp.append(ypp)
    print Files_Sign[i]
    i=i+3

################# pareggio tempi###############
aa=stats.mode(IYP)
for i in range (0, len(IYP)):
    delta_t=int(aa[0]-IYP[i])
    print (i,delta_t)
    if delta_t>0:
        Yp[i]= list([0] * delta_t)+list(Yp[i])
        Yp[i]=Yp[i][0:-delta_t]
        Ypp[i]= list([0] * delta_t)+list(Ypp[i])
        Ypp[i]=Ypp[i][0:-delta_t]
    elif delta_t<0:
        Yp[i]=Yp[i][abs(delta_t)-1:-1]
        Yp[i]= list(Yp[i])+list([0] *abs(delta_t))
        Ypp[i]=Ypp[i][abs(delta_t)-1:-1]
        Ypp[i]= list(Ypp[i])+list([0] *abs(delta_t))
    else:
        pass


IYN=[]    
i=0    
while i<=(len(Files_Sign)-1):
    yn=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=A)
    ynn=yn
    yn=butter_bandpass_filter(yn, lowpass, highpass, fs, order=Or)
    yn=yn[:int(len(yn)/cut)]
    yn=yn-(one*(mean(yn)))
    Yn.append(yn/max(abs(yn)))
    
       
    yn4=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=4)
    yn4=yn4/max(abs(yn4))
    iyn=pos(abs(yn4),0.02)
    IYN.append(iyn)
    #Ynf.append(butter_bandpass_filter(ynn, 70, 80, fs, order=Or))

    xnn=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=B)
    ynn=(ynn**2 + xnn**2)**0.5 *sign(ynn+xnn)
    ynn=ynn[:int(len(ynn)/cut)]
    ynn=ynn-(one*(mean(ynn)))
    ynn=butter_bandpass_filter(ynn, lowpass, highpass, fs, order=Or)
    Ynn.append(ynn)
    
    
    print Files_Sign[i]
    i=i+3

################# pareggio tempi###############
aa=stats.mode(IYN)
for i in range (0, len(IYN)):
    delta_t=int(aa[0]-IYN[i])
    print (i,delta_t)
    if delta_t>0:
        Yn[i]= list([0] * delta_t)+list(Yn[i])
        Yn[i]=Yn[i][0:-delta_t]
        Ynn[i]= list([0] * delta_t)+list(Ynn[i])
        Ynn[i]=Ynn[i][0:-delta_t]
    elif delta_t<0:
        Yn[i]=Yn[i][abs(delta_t)-1:-1]
        Yn[i]= list(Yn[i])+list([0] *abs(delta_t))
        Ynn[i]=Ynn[i][abs(delta_t)-1:-1]
        Ynn[i]= list(Ynn[i])+list([0] *abs(delta_t))
    else:
        pass    
    

i=1   
while i<=(len(Files_Sign)-1):
    xp=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=B)
    xp=butter_bandpass_filter(xp, lowpass, highpass, fs, order=Or)
    xp=xp[:int(len(xp)/cut)]
    xp=xp-(one*(mean(xp))) 
    Xp.append(xp/max(abs(xp)))
    #Xpf.append(butter_bandpass_filter(xpp, 70, 80, fs, order=Or))
    print Files_Sign[i]
    i=i+3

################# pareggio tempi###############
aa=stats.mode(IYP)
for i in range (0, len(IYP)):
    delta_t=int(aa[0]-IYP[i])
    if delta_t>0:
        Xp[i]= list([0] * delta_t)+list(Xp[i])
        Xp[i]=Xp[i][0:-delta_t]
        
    elif delta_t<0:
        Xp[i]=Xp[i][abs(delta_t)-1:-1]
        Xp[i]= list(Xp[i])+list([0] *abs(delta_t))
        
    else:
        pass


    
i=0    
while i<=(len(Files_Sign)-1):
    xn=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=B)
    xn=butter_bandpass_filter(xn, lowpass, highpass, fs, order=Or)
    xn=xn[:int(len(xn)/cut)]
    xn=xn-(one*(mean(xn)))
    Xn.append(xn/max(abs(xn)))#percentile(abs(z),84))
    #Xnf.append(butter_bandpass_filter(xnn, 70, 80, fs, order=Or))#percentile(abs(z),84))
    print Files_Sign[i]
    i=i+3
    
################# pareggio tempi###############
aa=stats.mode(IYN)
for i in range (0, len(IYN)):
    delta_t=int(aa[0]-IYN[i])
    if delta_t>0:
        Xn[i]= list([0] * delta_t)+list(Xn[i])
        Xn[i]=Xn[i][0:-delta_t]
    elif delta_t<0:
        Xn[i]=Xn[i][abs(delta_t)-1:-1]
        Xn[i]= list(Xn[i])+list([0] *abs(delta_t))
    else:
        pass    
  





for i in change:
    print i
    Yp[i]=Xp[i]
    Yn[i]=Xn[i]


Yp=array(Yp)
Yn=array(Yn)

Xp=array(Xp)
Xn=array(Xn)
Z=array(Z)

Ypf=array(Ypf)
Ynf=array(Ynf)

Xpf=array(Xpf)
Xnf=array(Xnf)

Ypp=array(Ypp)
Ynn=array(Ynn)

Xpp=array(Xpp)
Xnn=array(Xnn)



time=(arange(0,len(Z[0]), 1)*DT)/1000


YYp=(Yp*1.)
YYn=(Yn*1.)


#YYp=butter_bandpass_filter(YYp, lowpass, highpass-50, fs, order=Or)
#YYn=butter_bandpass_filter(YYn, lowpass, highpass-50, fs, order=Or)

grad_time= (((len(time)*dz)/Vmed**0.98)-(delay_time/DT))/(len(time))

one=ones(len(time))
timex=time+(0.001*one)
YYYp=[]
W=[]

YpF=[]
YnF=[]

        
for i in range (0, len(Yp)):
    w=((exp(-(timex-(delay_time+(grad_time*i*dz))/timex)**2)))
    w=w*(1-((exp(-(timex-(delay_time+(grad_time*i*dz))/timex)**6))))
    w=w/amax(w)
    print (i**0.9*grad_time) + delay_time
    
    YYn[i]=YYn[i]*w
    YYp[i]=YYp[i]*w
    YYp[i]=YYp[i]-YYn[i]
    YYn[i]=YYn[i]-YYp[i]
    YYp[i]=YYp[i]/amax(abs(YYp[i]))
    YYn[i]=YYn[i]/amax(abs(YYn[i]))
    
    #YYYp.append(YYp[i]-YYn[i])
    W.append(w)
    
   
   

#YYp=YYp-YYn
#YYn=YYn-YYp

time=(arange(0,len(Z[0]), 1)*DT)/1000

################################### Z  ##########################################

XX=[]

TT=[]

fig, az= plt.subplots(1,(int(len(Z))), sharey=True, facecolor='w', edgecolor='k',figsize=(35,10))

fig.subplots_adjust(left=.03, bottom= .05, right=.98, top=.95, hspace = .001, wspace=.001)

az = az.ravel()


for i in range (0, len(Z)):
    XX.append(0.)
    TT.append(0.)
    
for i in range(0, len(Z)):
        az[i].plot(Z[i][:(len(z))],time[:(len(z))], label='P-Waves')
        az[i].plot(XX[i],TT[i],marker='o',label='Picking')
        az[i].set_title("%dm"%((i+1)*dz),fontsize=9)
        az[i].yaxis.set_major_locator(ticker.MultipleLocator(0.01/cut))



pickz=[]
position=[]
prof=[]


def onclick(event):
    global iz, ix, iy
    iz,ix,iy = event.inaxes.get_title, event.xdata, event.ydata
    button=event.button
    x=event.xdata
    y=event.ydata
    A=("%s"%(iz))
    print A[-6:-2]
    
    if len(position)>=1:
        if "%s"%(A[-7:-2])in (position[::-1]):
            ii=1+position.index('%s'%(A[-7:-2]))
            prof.append(ii)
        else:
            position.append(A[-7:-2])
            ii=prof[-1]+1
            prof.append(ii)
    else:
        ii=1
        prof.append(ii)
        position.append(A[-7:-2])
        
    
    TT[(int(prof[-1]))-1]=y
    #fig.canvas.draw()
    print 'z(m) = %s, x(s) = %f'%(
       prof[-1], iy)
    
    O=(int(prof[-1]))-1
    az[O].plot(XX[O],TT[O],marker='o',color="r" )
    fig.canvas.draw()
                           
    global coords
    coords = [prof[-1], iy]
    
    pickz.append(coords)
    return coords

    # Refresh the plot
    
    #figure.canvas.draw()

cid = fig.canvas.mpl_connect('button_release_event', onclick)

    
plt.subplots_adjust(wspace=.0)
fig.canvas.set_window_title("Picking on P-Wave")

plt.show()
fig.savefig('OUTPUT\Picking on P-Waves.png')


pickz = sorted(pickz, key=lambda a_entry: a_entry[0])
savetxt('OUTPUT/pick_Z.asc',transpose(transpose(pickz)), fmt="%1.6f")
   

TTz=TT
profz=prof
################################## Y  ##########################################


XX=[]

TT=[]




fig1, ay = plt.subplots(1,(int(len(Yp))), sharey=True, facecolor='w', edgecolor='k')
fig1.subplots_adjust(left=.03, bottom= .05, right=.98, top=.95, hspace = .001, wspace=.001)


ay = ay.ravel()


for i in range (0, len(Z)):
    XX.append(0.)
    TT.append(0.)

for i in range(0, len(Yp)):
    ay[i].plot(Yp[i],time )
    ay[i].plot(Yn[i],time )
    ay[i].plot(XX[i],TT[i],marker='o' )
    #ay[i].plot(Yp[i]-Yn[i],time,linewidth=2.0, color= "red" )
    #ay[i].plot((Yp[i]-Yn[i]),time,'--', color= "grey" )
    ay[i].yaxis.set_major_locator(ticker.MultipleLocator(0.01/cut))
    if (i+1)*dz in (profz):
        t=profz.index((i+1)*dz)
        print t
        ay[i].plot(XX[t]*0.,TTz[t], marker='x', color= "black" )
    else:
        ay[i].plot(0., 0., marker='x', color= "black" )
    ay[i].set_title("%dm"%((i+1)*dz),fontsize=9)
    if i in (Zoom-ones(len(Zoom))):
        ay[i].set_xlim([-X_Zoom, X_Zoom])
    else:
        ay[i].set_xlim([-1, 1])

picky=[]
position=[]
prof=[]

def onclick(event):
    global iz, ix, iy
    iz,ix,iy = event.inaxes.get_title, event.xdata, event.ydata
    button=event.button
    x=event.xdata
    y=event.ydata
    A=("%s"%(iz))

    
    if len(position)>=1:
        if "%s"%(A[-7:-2])in (position[::-1]):
            ii=1+position.index('%s'%(A[-7:-2]))
            prof.append(ii)
        else:
            position.append(A[-7:-2])
            ii=prof[-1]+1
            prof.append(ii)
    else:
        ii=1
        prof.append(ii)
        position.append(A[-7:-2])
        
    
    TT[(int(prof[-1]))-1]=y
    #fig.canvas.draw()
    print 'z(m) = %s, x(s) = %f'%(
       prof[-1], iy)
    
    O=(int(prof[-1]))-1
    ay[O].plot(XX[O],TT[O],marker='o',color="r" )
    fig1.canvas.draw()
                           
    global coords
    coords = [prof[-1], iy]
    
    picky.append(coords)
    return coords
    
    # Refresh the plot
    
    #figure.canvas.draw()

cid = fig1.canvas.mpl_connect('button_release_event', onclick)

    
plt.subplots_adjust(wspace=.0)



plt.show()

fig1.savefig('OUTPUT\Picking on S-Waves.png')
picky = sorted(picky, key=lambda a_entry: a_entry[0])
savetxt('OUTPUT/pick_Y.asc',transpose(transpose(picky)), fmt="%1.6f")

################################## X  ##########################################

XX=[]

TT=[]

fig, ax = plt.subplots(1,(int(len(Yp))), sharey=True, facecolor='w', edgecolor='k')
fig.subplots_adjust(left=.03, bottom= .05, right=.98, top=.95, hspace = .001, wspace=.001)

ax = ax.ravel()

for i in range (0, len(Xp)):
    XX.append(0.)
    TT.append(0.)

for i in range(0, len(Xp)):
    ax[i].plot(XX[i],TT[i],marker='o' )
    ax[i].plot(YYp[i],time,linewidth=1.0, color= "c" )
    ax[i].plot(YYn[i],time,linewidth=1.0, color= "m" )
    ax[i].plot(W[i],time,':',linewidth=1.5, color= "black" )
    #ax[i].plot(Z[i],time,'--',linewidth=1.0, color= "grey" )
    ax[i].yaxis.set_major_locator(ticker.MultipleLocator(0.01/cut))
    if (i+1)*dz in (profz):
        t=profz.index((i+1)*dz)
        print t
        ax[i].plot(XX[t]*0.,TTz[t], marker='x', color= "black" )
    else:
        ax[i].plot(0., 0., marker='x', color= "black" )
    ax[i].set_title("%dm"%((i+1)*dz),fontsize=9)
    if i in (Zoom-ones(len(Zoom))):
        ax[i].set_xlim([-X_Zoom, X_Zoom])
    else:
        ax[i].set_xlim([-1, 1])
    
    
    ax[i].set_title("%dm"%((i+1)*dz),fontsize=9)


pickx=[]
positionx=[]
profx=[]

def onclick(event):
    global iz, ix, iy
    iz,ix,iy = event.inaxes.get_title, event.xdata, event.ydata
    button=event.button
    x=event.xdata
    y=event.ydata
    AA=("%s"%(iz))
    
    if len(positionx)>=1:
        if "%s"%(AA[-7:-2])in (positionx[::-1]):
            ii=1+positionx.index('%s'%(AA[-7:-2]))
            profx.append(ii)
        else:
            positionx.append(AA[-7:-2])
            ii=profx[-1]+1
            profx.append(ii)
    else:
        ii=1
        profx.append(ii)
        positionx.append(AA[-7:-2])

    TT[(int(profx[-1]))-1]=y
    fig.canvas.draw()
    print 'z(m) = %s, x(s) = %f'%(
       profx[-1], iy)
    
    O=(int(profx[-1]))-1
    ax[O].plot(XX[O],TT[O],marker='o',color="r" )
    fig.canvas.draw()
                           
    global coords
    coords = [profx[-1], iy]
    
    pickx.append(coords)
    return coords

    # Refresh the plot
    
    #figure.canvas.draw()

cid = fig.canvas.mpl_connect('button_release_event', onclick)

    
plt.subplots_adjust(wspace=.0)
fig.canvas.set_window_title("Picking on Treated S-Waves")

plt.show()

 

fig.savefig('OUTPUT\Picking on Treated S-Waves.png')

pickx = sorted(pickx, key=lambda a_entry: a_entry[0])
savetxt('OUTPUT/pick_X.asc',transpose(transpose(pickx)), fmt="%1.6f")



fig, ax = plt.subplots(1,int(len(YpF)), sharey=True, facecolor='w', edgecolor='k')
fig.subplots_adjust(left=.03, bottom= .05, right=.98, top=.95, hspace = .001, wspace=.001)

ax = ax.ravel()


for i in range(0, len(YpF)):
    #ax[i].plot(XX[i-1],TT[i-1],marker='o' )
    ax[i].plot(YpF[i],time,linewidth=1.5, color= "c" )
    #ax[i].plot(YnF[i],time,linewidth=1.0, color= "m" )
    #ax[i].plot(WW[i],time,':',linewidth=1.5, color= "black" )
    #ax[i].plot(Z[i],time,'--',linewidth=1.0, color= "grey" )
    ax[i].yaxis.set_major_locator(ticker.MultipleLocator(0.01/cut))
  

plt.show()



