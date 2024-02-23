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
    return popt 


def findfit(i1):
    cha=i1*1.0/100
    data=np.loadtxt("0-"+str('%.2f' % cha)+"-BAND.dat", delimiter=None)
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
    tfit=fit(Edata)
    hw=[]
    for i in range(300):
        if i<100:
            kx = np.pi*(99-i)/99
            ky = 0
        elif i<200:
            kx = np.pi*(i-100)/99; ky = kx
        else:
            kx = np.pi; ky = np.pi*(299-i)/99
        hw.append(eny4(tfit,kx,ky))
    return np.append(cha,tfit)
    

    
def mytry2():
    out=list()
    for i1 in range(80, 120):
        cha=i1*1.0/100
        data=np.loadtxt("0-"+str('%.2f' % cha)+"-BAND.dat", delimiter=None)
        press = np.zeros(4, dtype=np.int)
        a=0
        for i in range(0, len(data)-2):
            if data[i, 1]*data[i+1, 1]<0:
                if i%100!=99:
                    press[a]=(i//300)
                    a=a+1
        out.append(press)
    np.savetxt("femi.data", out, fmt="%.2f")

def mytry1():
    a1 = np.loadtxt('1.txt'); a2 = np.loadtxt('2.txt');a3 = np.loadtxt('3.txt');
    Edata = np.append(a1[100:200], a2[100:200])
    tfit=fit(Edata)
    hw=[]
    for i in range(300):
        if i<100:
            kx = np.pi*(99-i)/99
            ky = 0
        elif i<200:
            kx = np.pi*(i-100)/99; ky = kx
        else:
            kx = np.pi; ky = np.pi*(299-i)/99
        hw.append(eny4(tfit,kx,ky))
    print(hw)
    x = np.linspace(1, 300, 300)
    plt.plot(x,hw)
    plt.scatter(x,a1)
    plt.scatter(x,a2)
    plt.scatter(x,a3)
    plt.axis([1,300,-6,6]) 
    plt.show()

def mytry3():
    out=list()
    for i1 in range(80, 121):
        press=findfit(i1)
        out.append(press)
    np.savetxt("t.data", out, fmt="%.2f")

def mytry(i1):
    cha=i1*1.0/100
    #data=np.loadtxt("/home/ljcj007/vasp/20190702/1/0-"+str('%.2f' % cha)+"-BAND.dat", delimiter=None)
    data=np.loadtxt("/home/ljcj007/vasp/20191014_layers/"+str(i1)+"/1/BAND.dat", delimiter=None)
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
    tfit=fit(Edata)
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
    print(tfit)
    x = np.linspace(1, 300, 300)
    x2 = np.linspace(1, 300, 60)
    real=np.reshape(data[:, 1], (300,-1), order='F')  
    for i1 in range(np.shape(real)[1]):
        if i1%2==1:
            real[:, i1]=real[::-1, i1]
    return x,real,x2,hw
    
def fig4():
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(15,5))
    ax1 = plt.subplot(141)
    ax2 = plt.subplot(142)
    ax3 = plt.subplot(143)
    ax4 = plt.subplot(144)
    
    x,real,x2,hw=mytry(120)
    plt.sca(ax1)
    plt.plot(x2, hw, '^')
    plt.plot(x,real)
    plt.axis([1,300,-4,4]) 
    plt.xticks([0, 100, 200, 300],['X','$\Gamma$','M','X'])
    plt.grid(axis="x")    
    plt.grid(axis="y")    
    plt.ylabel('$E(meV)$')
    plt.title('(a)  $a/a_0$=1.2', loc ='left',fontsize='x-large', fontweight='heavy')
    
    x,real,x1,hw=mytry(110)
    plt.sca(ax2)
    plt.plot(x2, hw, '^')
    plt.plot(x,real)
    plt.axis([1,300,-4,4]) 
    plt.xticks([0, 100, 200, 300],['X','$\Gamma$','M','X'])
    plt.grid(axis="x")    
    plt.grid(axis="y")    
    plt.title('(b)  $a/a_0$=1.1', loc ='left',fontsize='x-large', fontweight='heavy')
    
    x,real,x1,hw=mytry(100)
    plt.sca(ax3)
    plt.plot(x2, hw, '^')
    plt.plot(x,real)
    plt.axis([1,300,-4,4]) 
    plt.xticks([0, 100, 200, 300],['X','$\Gamma$','M','X'])
    plt.grid(axis="x")    
    plt.grid(axis="y")    
    plt.title('(c)  $a/a_0$=1.0', loc ='left',fontsize='x-large', fontweight='heavy')
    
    
    x,real,x1,hw=mytry(90)
    plt.sca(ax4)
    plt.plot(x2, hw, '^')
    plt.plot(x,real)
    plt.axis([1,300,-4,4]) 
    plt.xticks([0, 100, 200, 300],['X','$\Gamma$','M','X'])
    plt.grid(axis="x")    
    plt.grid(axis="y")  
    plt.title('(d)  $a/a_0$=0.9', loc ='left',fontsize='x-large', fontweight='heavy')
    
    plt.show()

def fig1():
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(5,5))
    x,real,x2,hw=mytry(9)
    plt.plot(x2, hw, '^')
    plt.plot(x,real)
    plt.axis([1,300,-4,4]) 
    plt.xticks([0, 100, 200, 300],['X','$\Gamma$','M','X'],fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.grid(axis="x")    
    plt.grid(axis="y")    
    plt.ylabel('$E(eV)$',fontsize='x-large')
    plt.subplots_adjust(top=0.95, bottom=0.05, right=0.95, left=0.15, hspace=0.05, wspace=0.05)
    plt.margins(0, 0)
    plt.show()
    
def fig():
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(5,5))
    x,real,x2,hw=mytry(9)
    plt.plot(x2, hw, '^')
    plt.plot(x,real)
    plt.axis([150,250,-6,2]) 
    plt.xticks([200], ['M'], fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.grid(axis="x")      
    plt.ylabel('$Energy(eV)$',fontsize='x-large')
    #plt.subplots_adjust(top=0.95, bottom=0.05, right=0.95, left=0.25, hspace=0.05, wspace=0.05)
    plt.margins(0, 0)
    plt.show()
    
 
mytry3()
