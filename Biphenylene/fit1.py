import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from cmath import exp
    
def eny4(tpara,kx,ky):
    x1=0.193
    y1=0.338
    x2=0.500
    y2=0.160
    t1 = tpara[0]; t2 = tpara[1]; t3 = tpara[2]; t4 = tpara[3]
    s=np.zeros((6,6), dtype=complex)
    s[0,1] = -1*t1*exp((0.+1.j)*(ky*(1-2*y1))); 
    s[1,2] = -1*t2*exp((0.+1.j)*(kx*(x2+x1-1)+ky*(y1-y2)));
    s[2,3] = -1*t2*exp((0.+1.j)*(kx*(x2+x1-1)+ky*(y2-y1)));
    s[3,4] = -1*t1*exp((0.+1.j)*(ky*(2*y1-1))); 
    s[4,5] = -1*t2*exp((0.+1.j)*(kx*(1-x2-x1)+ky*(y1-y2)));
    s[5,0] = -1*t2*exp((0.+1.j)*(kx*(1-x2-x1)+ky*(y2-y1)));
    s[0,4] = -1*t1*exp((0.+1.j)*(kx*(2*x1))); 
    s[1,3] = -1*t1*exp((0.+1.j)*(kx*(2*x1))); 
    s[2,5] = -1*t3*exp((0.+1.j)*(ky*(2*y2))); 
    #s[0,2] = -1*t3*exp((0.+1.j)*(kx*(x1-x2)+ky*(1-y1-y2)));
    #s[1,3] = -1*t3*exp((0.+1.j)*(kx*(2*x1-1))); 
    #s[2,4] = -1*t3*exp((0.+1.j)*(kx*(x1-x2)+ky*(y1+y2-1)));
    #s[3,5] = -1*t3*exp((0.+1.j)*(kx*(x2-x1)+ky*(y1+y2-1)));
    #s[4,0] = -1*t3*exp((0.+1.j)*(kx*(1-2*x1))); 
    #s[5,1] = -1*t3*exp((0.+1.j)*(kx*(x2-x1)+ky*(1-y1-y2)));
    #s[0,3] = -1*t4*exp((0.+1.j)*(kx*(2*x1)+ky*(1-2*y1))); 
    #s[1,4] = -1*t4*exp((0.+1.j)*(kx*(2*x1)+ky*(2*y1-1))); 
    s=np.mat(s)
    s=s+s.H
    E,_=np.linalg.eig(s) 
    t=E.real
    return sorted(t)
    
    
def engy(K, t1, t2, t3, t4):
    nk=len(K)
    an=np.zeros(nk)
    tpara=[t1, t2, t3,t4]
    for i in range(0, nk-1):
        if i<50:
            kx = np.pi*(i)/49
            ky = 0
            t=eny4(tpara,kx,ky)
            an[i]=t[2]
        elif i<100:
            kx = np.pi; ky = np.pi*(i-50)/49
            t=eny4(tpara,kx,ky)
            an[i]=t[2]
        elif i<150:
            kx = np.pi*(149-i)/49; ky = kx
            t=eny4(tpara,kx,ky)
            an[i]=t[2]
        elif i<200:
            ky = np.pi*(i-150)/49; kx = 0
            t=eny4(tpara,kx,ky)
            an[i]=t[2]
        elif i<250:
            kx = np.pi*(i-200)/49; ky = np.pi
            t=eny4(tpara,kx,ky)
            an[i]=t[2]
        elif i<300:
            kx = np.pi*(i-250)/49
            ky = 0
            t=eny4(tpara,kx,ky)
            an[i]=t[3]
        elif i<350:
            kx = np.pi; ky = np.pi*(i-300)/49
            t=eny4(tpara,kx,ky)
            an[i]=t[3]
        elif i<400:
            kx = np.pi*(299-i)/49; ky = kx
            t=eny4(tpara,kx,ky)
            an[i]=t[3]
        elif i<450:
            ky = np.pi*(i-400)/49; kx = 0
            t=eny4(tpara,kx,ky)
            an[i]=t[3]
        elif i<500:
            kx = np.pi*(i-250)/49; ky = np.pi
            t=eny4(tpara,kx,ky)
            an[i]=t[3]
    return an 


def findfit(data):    
    real=np.reshape(data[:, 1], (250,-1), order='F')  
    for i1 in range(np.shape(real)[1]):
        if i1%2==1:
            real[:, i1]=real[::-1, i1]
    Edata = np.append(real[:, 11], real[:, 12])
    kxG_M = np.linspace(0,500,500)
    popt, pcov = curve_fit(engy, kxG_M, Edata, maxfev=50000, p0=(2.7,1.7,2.7,1.8))
    perr = np.sqrt(np.diag(pcov))
    print(popt)
    print(perr)
    return popt,perr

def engy1(K, t1, t2, t3, t4):
    nk=len(K)
    an=np.zeros(nk)
    tpara=[t1, t2, t3,t4]
    for i in range(0, nk-1):
        if i<50:
            
            kx = np.pi*(49-i)/49; ky = kx
            t=eny4(tpara,kx,ky)
            an[i]=t[2]
        elif i<100:
            kx = np.pi*(99-i)/49; ky = kx
            t=eny4(tpara,kx,ky)
            an[i]=t[3]
    return an 
def findfit1(data):    
    real=np.reshape(data[:, 1], (250,-1), order='F')  
    for i1 in range(np.shape(real)[1]):
        if i1%2==1:
            real[:, i1]=real[::-1, i1]
    Edata = np.append(real[100:150, 11], real[100:150, 12])
    kxG_M = np.linspace(0,100,100)
    popt, pcov = curve_fit(engy1, kxG_M, Edata, maxfev=50000, p0=(2.7,1.7,2.7,1.8))
    perr = np.sqrt(np.diag(pcov))
    print(popt)
    print(perr)
    return popt,perr

def energy_bond(data):
    tfit,perr=findfit(data)
    #tfit=[2.9,2.8,2.8,2.8]
    hw=[]
    for i in range(250):
        if i<50:
            kx = np.pi*(i)/49
            ky = 0
        elif i<100:
            kx = np.pi; ky = np.pi*(i-50)/49
        elif i<150:
            kx = np.pi*(149-i)/49; ky = kx
        elif i<200:
            ky = np.pi*(i-150)/49; kx = 0
        else:
            kx = np.pi*(i-200)/49; ky = np.pi
        hw.append(eny4(tfit,kx,ky))
    x2 = np.linspace(1, 250, 250)   
    return x2,hw


def fig(i):
    #data=np.loadtxt("./wannier90_band.dat")
    #x2=data[0:497,0]/4.1632844*250
    #hw=np.reshape(data[:, 1], (497,-1), order='F')+2.849335
    #x2 = np.linspace(1, 250, 497)
    data=np.loadtxt("./BAND-%.2d.dat" % i, delimiter=None)
    x2,hw=energy_bond(data)
    data=np.loadtxt("./BAND-%.2d.dat" % i, delimiter=None)
    x = np.linspace(1, 250, 250)
    real=np.reshape(data[:, 1], (250,-1), order='F')  
    for i1 in range(np.shape(real)[1]):
        if i1%2==1:
            real[:, i1]=real[::-1, i1]
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(3,4))
    plt.plot(x2, hw, linestyle='--')
    plt.plot(x,real)
    plt.axis([1,250,-5,5]) 
    plt.xticks([0, 50, 100, 150, 200,250],['$\Gamma$','K','M','$\Gamma$','K`','M'],fontsize='x-large')
    plt.yticks([-4, -2, 0, 2, 4],fontsize='x-large')
    plt.grid(axis="x")    
    plt.grid(axis="y")    
    plt.margins(0, 0)
    plt.subplots_adjust(top=0.98, bottom=0.10, right=0.98, left=0.12, hspace=0.05, wspace=0.05)
    plt.savefig('./'+str('%.2d' % i)+'.eps',dpi=600,format='eps')
    plt.show()

def energy_bond1():
    tfit=[2.9,2.9,2.8,2.8]
    hw=[]
    for i in range(500):
        if i<100:
            kx = np.pi*(i)/99
            ky = 0
        elif i<200:
            kx = np.pi; ky = np.pi*(i-100)/99
        elif i<300:
            kx = np.pi*(299-i)/99; ky = kx
        elif i<400:
            ky = np.pi*(i-300)/99; kx = 0
        else:
            kx = np.pi*(i-400)/99; ky = np.pi
        hw.append(eny4(tfit,kx,ky)[3])
    x2 = np.linspace(1, 500, 500)  
    hw=np.reshape(hw, (500,-1), order='F')  
    print(np.shape(hw))
    plt.figure(figsize=(3,4))
    plt.plot(x2[:], hw)
    plt.axis([1,500,-5,5]) 
    plt.xticks([0, 100, 200, 300, 400,500],['$\Gamma$','K','M','$\Gamma$','K`','M'],fontsize='x-large')
    plt.yticks([-4, -2, 0, 2, 4],fontsize='x-large')
    

#energy_bond1()
fig(12)
