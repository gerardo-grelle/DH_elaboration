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
cut=1 # cut the trace windows e.g. 2 is equal to half

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

####################
def moving_average_fo(a, n) :
    ret = cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n




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
lowpass=5
highpass=500
highpassp=500
fs=10000
Or=1

dx=DT/1000

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
y_freq_shot=[]
while i<=(len(Files_Sign)-1):
    
    yp4=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=4)
    yp4_A=yp4
    
    #yp4_A=(yp4[1:-1]-yp4[0:-2])/(DT*100000)
    #yp4=yp4/max(abs(yp4))
    iyp=pos(abs(yp4),0.02)
    IYP.append(iyp)
    nn= len(yp4)
    fx= yp4
    Fk= rfft(fx)
    nu= (rfftfreq(nn,dx))
    INPUT=abs(Fk*dx)
    My=max(INPUT)
    Mx=nu[INPUT.argmax()]
    print Mx
    yp4_A=butter_bandpass_filter(yp4_A, lowpass, highpass, fs, order=Or)
   

    
    yp4_A=yp4_A[:int(len(yp4_A)/cut)]
    Ypff.append(yp4_A)
    #Ypf.append(butter_bandpass_filter(ypp, 70, 80, fs, order=Or))


   
    #################################################
   
    ypp=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=A)
    xpp=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=B)
    #ypp=(ypp**2 + xpp**2)**0.5 *sign(cos(arctan(xpp/(ypp+0.0001))*sin(arctan(xpp/(ypp+0.0001)))))
    #ypp=(ypp**2 + xpp**2)**0.5 *-sign(ypp)
    #ypp=butter_bandpass_filter(ypp, lowpass, highpass, fs, order=Or) 
    ypp=ypp[:int(len(ypp)/cut)]
    #ypp=(ypp[1:-1]-ypp[0:-2])/(DT*100000)
    ypp=ypp-(one*(mean(ypp)))
    #ypp=butter_bandpass_filter(ypp, lowpass, highpass, fs, order=Or)
    Ypp.append(ypp/max(abs(yp4_A)))
    print Files_Sign[i]
    i=i+3
    
    ##### f0 da shot in y ########################
    ffyy= fft(yp4)
    fftyy = fftfreq(len(ffyy), DT/1000)[:len(ffyy)//2]
    ffyy=2.0/len(ffyy)*abs(ffyy[0:len(ffyy)//2])

    ffyy=moving_average_fo(ffyy,n=2)
    
    i_maxy=fftyy[argmax(ffyy)]
    y_freq_shot.append(i_maxy)

    plt.plot(fftyy[0:1000], ffyy[0:1000])
    #plt.plot(fftxx[1:-300], ffxx[1:-299]/max(ffxx[1:-299]))
    
    plt.xscale('log')
    #plt.yscale('log')
plt.show()
    






################# pareggio tempi###############
aa=stats.mode(IYP)
for i in range (0, len(IYP)):
    delta_t=int(aa[0]-IYP[i])
    print (i,delta_t)
    if delta_t>0:
        Ypp[i]= list([0] * delta_t)+list(Ypp[i])
        Ypp[i]=Ypp[i][0:-delta_t]
    elif delta_t<0:
        Ypp[i]=Ypp[i][abs(delta_t)-1:-1]
        Ypp[i]= list(Ypp[i])+list([0]*abs(delta_t))
    else:
        pass


IYN=[]    
i=0

x_freq_shot=[]
while i<=(len(Files_Sign)-1):
    
    
       
    yn4=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=4)
    yn4_A=yn4
    yp4_A=butter_bandpass_filter(yp4_A, lowpass, highpass, fs, order=Or)
    yn4_A=(yn4[1:-1]-yn4[0:-2])/(DT*100000)
    #yn4=yn4/max(abs(yn4))
    iyn=pos(abs(yn4),0.02)
    IYN.append(iyn)

    yn4_A=butter_bandpass_filter(yn4_A, lowpass, highpass, fs, order=Or)
    yn4_A=yp4_A[:int(len(yn4_A)/cut)]
    Ynff.append(yn4_A)

    
    #################################################

    
    #Ynf.append(butter_bandpass_filter(ynn, 70, 80, fs, order=Or))
    ynn=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=A)
    xnn=genfromtxt("%s"%(Files_Sign[i]),skip_header=1, delimiter=",", usecols=B)
    #ynn=(ynn**2 + xnn**2)**0.5 *sign(cos(arctan(xnn/(ynn+0.0001))*sin(arctan(xnn/(ynn+0.0001)))))
    #ynn=(ynn**2 + xnn**2)**0.5 *-sign(ynn)
    #ynn=butter_bandpass_filter(ynn, lowpass, highpass, fs, order=Or)
    #ynn=(ynn**2 + xnn**2)**0.5# *sign(ynn)
    ynn=ynn[:int(len(ynn)/cut)]
    #ynn=(ynn[1:-1]-ynn[0:-2])/(DT*100000)
    ynn=ynn-(one*(mean(ynn)))
    #ynn=butter_bandpass_filter(ynn, lowpass, highpass, fs, order=Or)
    Ynn.append(ynn/max(abs(yn4_A)))
   
    
    
    print Files_Sign[i]
    i=i+3

    ##### f0 da shot in x ########################
    ffxx= fft(yn4)
    fftxx = fftfreq(len(ffxx), DT/1000)[:len(ffxx)//2]
    ffxx=2.0/len(ffxx)*abs(ffxx[0:len(ffxx)//2])
    
    ffxx=moving_average_fo(ffxx,n=2)
    
    i_maxx=fftxx[argmax(ffxx)]
    x_freq_shot.append(i_maxx)
    

    plt.plot(fftxx[0:1000], ffxx[0:1000])
    #plt.plot(fftxx[1:-300], ffxx[1:-299]/max(ffxx[1:-299]))
    
    plt.xscale('log')
    #plt.yscale('log')
plt.show()


    

################# pareggio tempi###############
aa=stats.mode(IYN)
for i in range (0, len(IYN)):
    delta_t=int(aa[0]-IYN[i])
    print (i,delta_t)
    if delta_t>0:
        Ynn[i]= list([0] * delta_t)+list(Ynn[i])
        Ynn[i]=Ynn[i][0:-delta_t]
    elif delta_t<0:
        Ynn[i]=Ynn[i][abs(delta_t)-1:-1]
        Ynn[i]= list(Ynn[i])+list([0] *abs(delta_t))
    else:
        pass    
    


mfy=mean(y_freq_shot)
mfx=mean(x_freq_shot)


for i in change:
    print i
    Yp[i]=Xp[i]
    Yn[i]=Xn[i]


Ypp=array(Ypp)
Ynn=array(Ynn)

Z=array(Z)


time=(arange(0,len(Z[0]), 1)*DT)/1000



VsS= genfromtxt("Vs_VP_coherent_interval.asc",skip_header=1, delimiter=",", usecols=2)#simple time
VsS=array(VsS)

#***************** import median first time by picking  **************************
ftime=genfromtxt("OUTPUT/Vs_median_picking.asc",skip_header=1, delimiter=",", usecols=1)
ftime=array(ftime)
#YYp=butter_bandpass_filter(YYp, lowpass, highpass-50, fs, order=Or)
#YYn=butter_bandpass_filter(YYn, lowpass, highpass-50, fs, order=Or)

grad_time= (((len(time)*dz)/Vmed**0.98)-(delay_time/DT))/(len(time))
one=ones(len(time))
timex=time+(0.001*one)
YYYp=[]
W=[]



YpF=[]
YnF=[]
WW=[]
Ycut=[]
Xcut=[]
Ntcut_i=[]
Ntcut_i.append(-10)
Li=[]
Lix=[]
W=[]
Lo=2.50
zo=(dz/cos(arctan(2.50/dz)))
dzi=[]
pI=[]#slowless lungo le interfacce
Rho=[]

############################## max spectral furior for each signal 


def filt_wind (fo):
    A0=(fo-1,fo+1)
    A1=(fo-3,fo+3)
    A2=(fo-5,fo+5)
    A3=(fo-7,fo+7)
    A4=(fo-9,fo+9)
    A5=(fo-11,fo+11)
    A6=(fo-13,fo+13)
    A7=(fo-15,fo+15)
    B1=(A1[0]-3,A1[0]+3)
    B2=(B1[1]-3,B1[1]+3)
    B3=(B2[1]-3,B2[1]+3)
    return (A0,A0,A0,A1,A1,A2,A2,A3,A4,B1,B2,B3)

def rep (A,n):
    AA=[]
    for i in range(0,n):
        AA.append(A)
    return AA

#define a moving average function
def moving_average(x,y,step_size=0.05,width=1):
    bin_centers  = arange(min(x),max(x)-0.5*step_size,step_size)+0.5*step_size
    bin_avg = zeros(len(bin_centers))

    #We're going to weight with a Gaussian function
    def gaussian(x,amp=1,mean=0,sigma=1):
        return amp*exp(-(x-mean)**2/(2*sigma**2))

    for index in range(0,len(bin_centers)):
        bin_center = bin_centers[index]
        weights = gaussian(x,mean=bin_center,sigma=width)
        bin_avg[index] = average(y,weights=weights)

    return (bin_centers,bin_avg)


for i in range (0, len(Ypp)-1):
    ########### refraction running of the seismic ray
    inc=arcsin(sin(arctan(2.50/((i+1)*dz)))*(VsS[i+1]/VsS[i]))
    dzi.append(dz/cos(inc))
    print dz/cos(inc)
    ########## Slowless lungo le interfacce
    pI.append((sin(arctan(2.50/((i+1)*dz)))*(VsS[i+1]/VsS[i]))/VsS[i])


    ########### Medium density
    Rho.append(4.4*VsS[i]**0.25)
T_SH=[]
for i in range (0, len(Ypp)-1):
    qi=1/VsS[i]*(1-(VsS[i]* pI[i])**2)**0.5
    qii=1/VsS[i]*(1-(VsS[i]* pI[i])**2)**0.5
    if i==0:
        t_SH_o=1.
        T_SH.append(t_SH_o)
    else:
        t_SH= (2*Rho[i-1]*(VsS[i-1]**2)*qi)/((Rho[i-1]*(VsS[i-1]**2)*qi)+(Rho[i]*(VsS[i]**2)*qii))
        T_SH.append(t_SH)

for i in range (0, len(Ypp)-1):
    print (i+1, VsS[i], T_SH[i]) 
  
   
omg=6.28*mfy
Timex=[]
for ii in range (1,6):
    Ycuto=[]
    Xcuto=[]
    for i in range (0, len(Ypp)):
        ############## REFRACTION RUNNING ############
       
        #lz=((i+1)*dz*cos(arctan(3.40/((i+1)*dz))))- (i*dz*cos(arctan(3.40/(i*dz))))
        #lii=((i+1)*dz*cos(arctan(3.40/((i+1)*dz))))
        li=(dz/cos(arctan(2.50/dz)))
        Li.append(li)
        lix=(dz*(i+1))/cos(arctan(2.5/(i+1)*dz))
        Lix.append(lix)
        ntcut_i=int(ftime[i]/(DT/1000))
        Ntcut_i.append(ntcut_i)
        ntcut_f= ntcut_i+int((0.5*ii/(mfy*DT))*1000)
        ycut=Ypp[i][ntcut_i:ntcut_f]
        xcut=Ynn[i][ntcut_i:ntcut_f]
        timexx=timex
        timexx=timexx[ntcut_i:ntcut_f]
        #fig, ax = plt.subplots()
        #ax.plot(timexx, ycut)


        
        timexx= (timexx-timexx[0])/(timexx[-1]-timexx[0])
        w=1-exp(-(timexx/0.1)**2)
        w=w*w[::-1]
        w=w/max(w)
        print (len(timexx),len(ycut))
       
        
       
        ycut_=[]
        xcut_=[]
        FILy= filt_wind(((mfy+1)) -(i*(dz+1))**0.0)
        FILx= filt_wind(((mfx+1)) -(i*(dz+1))**0.0)
        
        yfin=ycut
        yfin= pad(yfin, (0, 5000), 'constant')
        
        yfin= fft(yfin)
        fo_last = fftfreq(len(yfin), DT/1000)[:len(yfin)//2]
        yfin=2.0/len(yfin)*abs(yfin[0:len(yfin)//2])

        yfin=moving_average_fo(yfin,n=10)
            
        fo_lasti=fo_last[argmax(yfin)]
        print fo_lasti
        print ((mfy+1) -(i*(dz+1))**0.0)
        
        
    



        
        for j in range (0,len(FILy)):
            ycut=butter_bandpass_filter(ycut, FILy[j][0], FILy[j][1], fs, order=Or)
            xcut=butter_bandpass_filter(xcut, FILx[j][0], FILx[j][1], fs, order=Or)
            ycut_.append(ycut*w) 
            xcut_.append(xcut*w)
        Ycuto.append(ycut_)  
        Xcuto.append(xcut_)    
        
    Timex.append(timexx)        
    Ycut.append(Ycuto)
    Xcut.append(Xcuto)



#######################################################################################
App=[]
Bpp=[]

Ann=[]
Bnn=[]

Ayp=[]
Byp=[]
Ayn=[]
Byn=[]

dampTy=[]
dampTx=[]
#print (amax(abs(Ypff[i]))/max(amax(Ypff,1)))
for ii in range (0,5):
    Ayp1=[]
    Byp1=[]
    Ayn1=[]
    Byn1=[]
    FF=[]
    for i in range (1, len(Ypp)-1):
        I1=VsS[i-1]/VsS[i]
        #I2=VsS[i]/VsS[i+1]
        #F1= 1/((cos(omg*(dz/2)/VsS[i-1]))**2 +((1/I1**2)*(sin(omg*(dz/2)/VsS[i-1])**2)))**0.5
        F2= T_SH[i]
        F=F2#*F2
        FF.append(F)
        #print F1
        #print"++++++++++++++++++++++"
        #print i,VsS[i],F
        Ayp2=[]
        Byp2=[]
        Ayn2=[]
        Byn2=[]

        
        for j in range (0, len(FILy)):
            Ap=mean((abs(Ycut[ii][i-1][j]))/(Lix[i+1]))
            Bp=mean((abs(Ycut[ii][i][j]))/(FF[i-1]*(Lix[i])))
            Ayp2.append(((Ycut[ii][i-1][j]))/(Lix[i+1]))
            Byp2.append(((Ycut[ii][i][j]))/(FF[i-1]*(Lix[i])))
            
            AA=(log((Ap/Bp)/1.0001))/(((Lix[i+1]-Lix[i])*(mfy*2*pi))/(VsS[i]))
            
            An=mean((abs(Xcut[ii][i-1][j]))/(Lix[i+1]))
            Bn=mean((abs(Xcut[ii][i][j]))/(FF[i-1]*(Lix[i])))
            Ayn2.append(((Xcut[ii][i-1][j]))/(Lix[i+1]))
            Byn2.append(((Xcut[ii][i][j]))/(FF[i-1]*(Lix[i])))
            BB=(log((An/Bn)/1.0001))/(((Lix[i+1]-Lix[i])*(mfx*2*pi))/(VsS[i]))
            
            
            if -0.20<=AA<=0.20:
                print i, AA, VsS[i-1],VsS[i],FF[i-1]
                dampTx.append(AA)
                dampTy.append(-i)
            if -0.20<=BB<=0.20:
                print i, BB, VsS[i-1],VsS[i],FF[i-1]
                dampTx.append(BB)
                dampTy.append(-i)

        Ayp1.append(Ayp2)
        Byp1.append(Byp2)
        Ayn1.append(Ayn2)
        Byn1.append(Byn2)

    Ayp.append(Ayp1)
    Byp.append(Byp1)
    Ayn.append(Ayn1)
    Byn.append(Byn1)


# Function to convert  
def list_String(s): 
    str1 = [] 
    # traverse in the string  
    for ele in s: 
        str1.append("%d"%(ele))  
    # return string  
    return str1 




Damp_all=[]
Med_damp=[]
P66damp=[]
P33damp=[]
depth=[]
for i in range (0,len (VsS)):
    DampT=[]
    for j in range (1,len(dampTx)):
        if dampTy[j]==-i*dz:
            DampT.append(dampTx[j])
        else:
            pass
    if DampT!=[]:
        Med_damp.append(median(DampT))
        P33damp.append(percentile(DampT, 37.5))
        P66damp.append(percentile(DampT, 62.5))
        depth.append(-i*dz)
        Damp_all.append(DampT)                

lab=list_String(depth)
aa=abs(array(depth))
plt.rcParams["figure.figsize"] = (8,12)
plt.boxplot(Damp_all, vert=False)
plt.yticks(aa, lab)
plt.gca().invert_yaxis()

plt.show()
MM_damp=[]
P33_damp=[]
P66_damp=[]
Depth=[]
for i in range(0, len(Med_damp)):
    if isnan(Med_damp[i])==True:
        pass
    else:
        MM_damp.append(Med_damp[i])
        P33_damp.append(P33damp[i])
        P66_damp.append(P66damp[i])
        Depth.append(depth[i])


dampTy=array(dampTy)+1.
Depth=array(Depth)+1.
Med_mob2=moving_average(Depth,MM_damp,step_size=0.5,width=1)
Med_mob3=moving_average(Depth,MM_damp,step_size=1,width=2)
Med_mob3_33=moving_average(Depth,P33_damp,step_size=1,width=2)
Med_mob3_66=moving_average(Depth,P66_damp,step_size=1,width=2)
plt.rcParams["figure.figsize"] = (8,12)
plt.plot(dampTx,dampTy,'.',color='gray')
plt.plot(MM_damp,Depth,color='black')
plt.plot(Med_mob2[1],Med_mob2[0], color='red',  linewidth=2)
plt.plot(Med_mob3[1],Med_mob3[0], color='black',  linewidth=3)
plt.plot(Med_mob3_33[1],Med_mob3[0], color='black',  linewidth=3,linestyle='--')
plt.plot(Med_mob3_66[1],Med_mob3[0], color='black',  linewidth=3, linestyle='--')

plt.fill_betweenx(Med_mob3[0], Med_mob3_33[1], Med_mob3_66[1], alpha=0.2)
plt.xlim(-0.25, 0.25)
plt.grid(which='major', color='gray', linewidth=0.5)
plt.grid(which='minor', color='gray', linewidth=0.5)
mindepth=-(int(round(max(abs(Depth)))/5.0)*5.0)-5.
plt.ylim(mindepth,0)
plt.yticks(arange(mindepth,0 , 5))
plt.xticks(arange(-0.20, 0.225, 0.025))
plt.xticks(rotation = 90)
plt.plot([0,0],[mindepth,0],color='gold',linewidth=5,linestyle='dashed')
plt.show()


print (len([1 for i in dampTx if i > 0]))
print (len([1 for i in dampTx if i < 0]))

j=0.
AA=[-1.,1.]
while min(AA)<0. and j<0.005:
    AA=Med_mob3[1]+(j)
    j=j+0.00005
    
print j

CC=[4, 14.,21.,25.,30.,36.5,41.]
P_z=[]
P_Damp=[]
i=0
Med_damp_P=[]
zz=abs(Med_mob3[0][::-1])
CC.insert(len(CC),zz[-1])
dd=Med_mob3[1][::-1]+(j)
P_z.append(zz[0])
for ii in range (0,len(zz)):
    Med_damp_P.append(dd[ii])
    if int(CC[i])== int(zz[ii]):
        P_z.append(-zz[ii])
        P_z.append(-zz[ii])
        P_Damp.append(mean(Med_damp_P))
        P_Damp.append(mean(Med_damp_P))
        print Med_damp_P
        Med_damp_P=[]
        i=i+1
        
    else:
        pass
    
        
       
            
sc=1.-j**0.5       
        
#j=0   




P_Damp=array(P_Damp)

##################### PLOT CORRECTED DAMPING PROFILE ##############################
plt.rcParams["figure.figsize"] = (8,12)
plt.plot(Med_mob2[1]+(j)*sc,Med_mob2[0], color='red',  linewidth=1.0,linestyle='--')
plt.plot((Med_mob3[1]+(j))*sc,Med_mob3[0], color='black',  linewidth=1.5)
plt.plot(P_Damp*sc,P_z[0:-1], color='blue',  linewidth=3.5)
#plt.plot(Med_mob3_33[1]+(j*abs(Med_mob3[0])),Med_mob3[0], color='black',  linewidth=3,linestyle='--')
#plt.plot(Med_mob3_66[1]+(j*abs(Med_mob3[0])),Med_mob3[0], color='black',  linewidth=3, linestyle='--')

plt.fill_betweenx(Med_mob3[0], (Med_mob3_33[1]+(j))*sc, (Med_mob3_66[1]+(j))*sc, alpha=0.1)
plt.grid(which='major', color='gray', linewidth=0.5)
plt.grid(which='minor', color='gray', linewidth=0.5)
mindepth=-(int(round(max(abs(Depth)))/5.0)*5.0)-5.
plt.ylim(mindepth,0)
plt.yticks(arange(mindepth,0 , 5))
plt.xticks(arange(0.0, 0.15, 0.01))
#plt.plot([0,0],[mindepth,0],color='yellow',linewidth=3)
plt.xlim(0.0, 0.1)
plt.show()





plt.rcParams["figure.figsize"] = (8,5)

for i in range (0, len(FILy)):
    fig, ax = plt.subplots()
    
    ax.plot(Timex[0], Ayp[0][22][i])
    ax.plot(Timex[0], Byp[0][22][i])

    ax.set(xlabel='Relative Period (k/n)', ylabel='Amplitude',
           title='k=1; freq.range => [%d, %d]'%(FILx[i][0],FILx[i][1]))
    ax.grid()


    plt.show()










print len(fo_last)
print fo_last

x=xxxxxxx
for i in range (0, len(VsS)):
    print (i, AA, VsS[i+1],FF[i])            


        
        
        
        #App.append(Ap)
        #Bpp.append(Bp)
        

        #An=mean(abs(Xcut[ii][i-1]))/(sum(dzi[i]))
        #Bn=mean(abs(Xcut[ii][i]))/F/(sum(dzi[i-1]))
        #Ann.append(An)
        #Bnn.append(Bn)
        
        
        #print (log((A*Li[i-1])/(B*Li[i]))*VsS[i])/(75*6.28)
        #AA= (log((Ap/Bp)/1.001))/((dzi[i]*60*2*pi)/(VsS[i]))
        #BB=(log((An/Bn)/1.001))/((dzi[i]*60*2*pi)/(VsS[i]))
        #if -0.20<=AA<=0.20:
            #print i, AA
        #if -0.20<=BB<=0.20:
            #print i, BB



fig, ax = plt.subplots()
ax.plot(timexx, Ycut[2][13])
ax.plot(timexx, Ycut[2][13])

ax.set(xlabel='time (s)', ylabel='voltage (mV)',
       title='About as simple as it gets, folks')
ax.grid()


plt.show()


    
   
fig, ax = plt.subplots(1,3, sharey=True, facecolor='w', edgecolor='k')
fig.subplots_adjust(left=.03, bottom= .05, right=.98, top=.95, hspace = .001, wspace=.001)

ax = ax.ravel()




for i in range(7, 10):
    #ax[i].plot(XX[i-1],TT[i-1],marker='o' )
    ax[i].plot(Ycut[i],timexx,linewidth=1.5, color= "c" )
    ax[i].plot(Ycut[i],timexx,linewidth=1.5, color= "m" )
    #ax[i].plot(WW[i],time,':',linewidth=1.5, color= "black" )
    #ax[i].plot(Z[i],time,'--',linewidth=1.0, color= "grey" )

  

plt.show()





