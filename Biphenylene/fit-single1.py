import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from cmath import exp

def matr(tpara,kx,ky):
    x1=0.193
    y1=0.338
    x2=0.500
    y2=0.160
    t1 = tpara[0]; t2 = tpara[1]; t3 = tpara[2]; t4 = tpara[3]
    s=np.zeros((6,6), dtype=np.complex)
    s[0,1] = -1*t1*exp((0.+1.j)*(ky*(1-2*y1))); 
    s[1,2] = -1*t2*exp((0.+1.j)*(kx*(x2+x1-1)+ky*(y1-y2)));
    s[2,3] = -1*t2*exp((0.+1.j)*(kx*(x2+x1-1)+ky*(y2-y1)));
    s[3,4] = -1*t1*exp((0.+1.j)*(ky*(2*y1-1))); 
    s[4,5] = -1*t2*exp((0.+1.j)*(kx*(1-x2-x1)+ky*(y1-y2)));
    s[5,0] = -1*t2*exp((0.+1.j)*(kx*(1-x2-x1)+ky*(y2-y1)));
    s[0,4] = -1*t3*exp((0.+1.j)*(kx*(2*x1))); 
    s[1,3] = -1*t3*exp((0.+1.j)*(kx*(2*x1))); 
    s[2,5] = -1*t4*exp((0.+1.j)*(ky*(2*y2))); 
    s=np.mat(s)
    s=s+s.H
    return s
    
def eny4(tpara,kx,ky):
    E,_=np.linalg.eig(matr(tpara,kx,ky)) 
    t=E.real
    return sorted(t)
    
    
def engy(K, t1, t2, t3, t4):

    if K<50:
        kx = np.pi*(K)/49
        ky = 0
    elif K<100:
        kx = np.pi; ky = np.pi*(K-50)/49
    elif K<150:
        kx = np.pi*(199-K)/49; ky = kx
    elif K<200:
        ky = np.pi*(K-150)/49; kx = 0
    else:
        kx = np.pi*(K-200)/49; ky = np.pi
        
    tpara=[t1, t2, t3,t4]
    t=eny4(tpara,kx,ky)
    
    return an 
    
def fit(Edata):
    kxG_M = np.linspace(0, 500)
    popt, pcov = curve_fit(engy, kxG_M, Edata, maxfev=50000, p0=(-2.91,-3.19,-2.82,-2.80))
    perr = np.sqrt(np.diag(pcov))
    print(popt)
    print(perr)
    return popt,perr


def findfit(data):    
    a1=data[:, 11]
    a2=data[:, 12]
    a2=a2[::-1]
    Edata = np.append(a1[:], a2[:])
    popt,perr=fit(Edata)
    return popt,perr
    
def get_t():
    out=list()
    for i1 in range(90, 121):
        cha=i1*1.0/100    
        data=np.loadtxt("../20190702/1/0-"+str('%.2f' % cha)+"-BAND.dat", delimiter=None)
        popt,perr=findfit(data)
        press=np.append(cha,popt)
        press=np.append(press,perr)
        out.append(press)
    np.savetxt("t.data", out, fmt="%.4f")

def energy_bond(data):
    tfit,perr=findfit(data)
    hw=[]
    for iii in range(60):
        i=iii*5
        if i<100:
            kx = np.pi*(99-i)/99
            ky = 0
        elif i<200:
            kx = np.pi*(i-100)/99; ky = kx
        else:
            kx = np.pi; ky = np.pi*(299-i)/99
        hw.append(eny4(tfit,kx,ky))
    x = np.linspace(1, 300, 300)
    x2 = np.linspace(1, 300, 60)
    real=np.reshape(data[:, 1], (300,-1), order='F')  
    for i1 in range(np.shape(real)[1]):
        if i1%2==1:
            real[:, i1]=real[::-1, i1]
    return x,real,x2,hw


def fig(i):
    cha=i*1.0/100    
    data=np.loadtxt("../20190702/1/0-"+str('%.2f' % cha)+"-BAND.dat", delimiter=None)
    x,real,x2,hw=energy_bond(data)
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(3,4))
    plt.plot(x2, hw, linestyle='--')
    plt.plot(x,real)
    plt.axis([1,300,-4.5,4.5]) 
    plt.xticks([0, 100, 200, 300],['X','$\Gamma$','M','X'],fontsize='x-large')
    plt.yticks([-4, -2, 0, 2, 4],fontsize='x-large')
    plt.grid(axis="x")    
    plt.grid(axis="y")    
    plt.margins(0, 0)
    plt.text(140, 3.2, r'$\epsilon_1$',fontsize = 'xx-large') 
    plt.text(170, 0.5, r'$\epsilon_2$',fontsize = 'xx-large') 
    plt.text(110, -0.8, r'$\epsilon_3$',fontsize = 'xx-large') 
    plt.text(140, -3.2, r'$\epsilon_4$',fontsize = 'xx-large') 
    plt.subplots_adjust(top=0.98, bottom=0.10, right=0.98, left=0.12, hspace=0.05, wspace=0.05)
    plt.savefig('./'+str('%.2d' % i)+'.eps',dpi=600,format='eps')
    plt.show()

def energy_bond1():
    tfit=[2.45,2.45,2.45,2.45]
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
        hw.append(eny4(tfit,kx,ky))
    x2 = np.linspace(1, 500, 500)    
    plt.figure(figsize=(3,4))
    plt.plot(x2, hw, linestyle='--')
    plt.axis([1,500,-15,15]) 
    plt.xticks([0, 100, 200, 300, 400,500],['$\Gamma$','K','M','$\Gamma$','K`','M'],fontsize='x-large')
    plt.yticks([-4, -2, 0, 2, 4],fontsize='x-large')
    

energy_bond1()
