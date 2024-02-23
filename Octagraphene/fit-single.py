import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import cmath

def matr(tpara,kx,ky):
    a=0.3
    t = tpara[0]; t1 = tpara[1]; t2 = tpara[2]
    s=np.zeros((4,4), dtype=np.complex)
    s[0,1] = -1*t*cmath.exp((0.+1.j)*a*(kx-ky)); 
    s[0,2] = -1*t1*cmath.exp((0.+1.j)*(1-2*a)*ky)-t2*cmath.exp(-(0.+1.j)*2*a*ky);
    s[0,3] = -1*t*cmath.exp((0.+1.j)*a*(-kx-ky));
    s[1,2] = -1*t*cmath.exp(-(0.+1.j)*a*(kx+ky));
    s[1,3] = -1*t1*cmath.exp((0.+1.j)*(1-2*a)*kx)-t2*cmath.exp(-(0.+1.j)*2*a*kx);
    s[2,3] = -1*t*cmath.exp((0.+1.j)*a*(-kx+ky));
    s=np.mat(s)
    return s+s.H
    
def eny4(tpara,kx,ky):
    E,_=np.linalg.eig(matr(tpara,kx,ky)) 
    t=E.real
    return sorted(t)
    
    
def engy(K, t, t1, t2):
    kxx = K
    kyy = K
    tpara=[t, t1, t2]
    nk=len(kxx)
    an=np.zeros(nk)
    for i in range(0, nk-1):
        kx = kxx[i]; ky = kyy[i]
        t=eny4(tpara,kx,ky)
        if i<100:
            an[i] = t[2]
        elif i<200:
            an[i] = t[1]
        else:
            an[i] = t[0]
    return an 
    
def fit(Edata):
    #kxX_G=np.linspace(np.pi, 0,100)
    #kyX_G=0*kxX_G
    kxG_M = np.linspace(0, np.pi, 100)
    #kyG_M=kxG_M
    #kxM_X = np.pi*np.ones(100)
    #kyM_X=np.linspace(np.pi, 0,100)
    K=np.append(kxG_M,kxG_M)
    popt, pcov = curve_fit(engy, K, Edata, maxfev=50000, p0=(-2.91,-3.19,-0.82))
    perr = np.sqrt(np.diag(pcov))
    print(popt)
    print(perr)
    return popt,perr


def findfit(data):
    press = np.zeros(4, dtype=np.int)
    a=0
    for i in range(0, len(data)-2):
        if data[i, 1]*data[i+1, 1]<0:
            if i%100!=99:
                press[a]=(i//300)
                a=a+1
    if press[0]==press[1]:
        a2=press[0]
        a1=press[2]
    else:
        a2=press[0]
        a1=press[1]
    a1=data[a1*300:a1*300+300, 1]
    a2=data[a2*300:a2*300+300, 1]
    a2=a2[::-1]
    Edata = np.append(a1[100:200], a2[100:200])
    popt,perr=fit(Edata)
    return popt,perr
    
def get_t():
    out=list()
    for i1 in range(90, 121):
        cha=i1*1.0/100    
        data=np.loadtxt("1/0-"+str('%.2f' % cha)+"-BAND.dat", delimiter=None)
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
    data=np.loadtxt("1/0-"+str('%.2f' % cha)+"-BAND.dat", delimiter=None)
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

fig(1.00)
