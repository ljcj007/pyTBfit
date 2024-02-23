import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import cmath    
def matr(tpara,kx,ky):
    t = tpara[0]; t1 = tpara[1]; t2 = tpara[2]; t3 = tpara[3]
    s=np.zeros((4,4), dtype=np.complex)
    s[0,1] = -1*t; 
    s[0,2] = -1*t1*cmath.exp((0.+1.j)*ky)-t2;
    s[0,3] = -1*t;
    s[1,2] = -1*t;
    s[1,3] = -1*t1*cmath.exp((0.+1.j)*kx)-t2;
    s[2,3] = -1*t;
    s=np.mat(s)
    s=s+s.H
    s1=-1*t3*np.identity(4)
    s0=0*np.identity(4)
    v1=np.vstack((np.vstack((s, s1)), np.vstack((s0, s0))))
    v2=np.vstack((np.vstack((s1, s)), np.vstack((s1, s0))))
    v3=np.vstack((np.vstack((s0, s1)), np.vstack((s, s1))))
    v4=np.vstack((np.vstack((s0, s0)), np.vstack((s1, s))))
    E,_=np.linalg.eig(np.hstack((np.hstack((v1, v2)), np.hstack((v3, v4)))))
    t=E.real
    return sorted(t)
    
def engy(K, t, t1, t2, t3):
    kxx = K
    kyy = K
    tpara=[t, t1, t2, t3]
    nk=len(kxx)
    an=np.zeros(nk)
    for i in range(0, nk-1):
        kx = kxx[i]; ky = kyy[i]
        t=matr(tpara,kx,ky)
        if i<100:
            an[i] = t[4]
        elif i<200:
            an[i] = t[5]
        elif i<300:
            an[i] = t[6]
        elif i<400:
            an[i] = t[7]
        elif i<500:
            an[i] = t[8]
        elif i<600:
            an[i] = t[9]
        elif i<700:
            an[i] = t[10]
        elif i<800:
            an[i] = t[11]
        else:
            an[i] = t[0]
    return an 
    
def fit(Edata):
    kxG_M = np.linspace(0, np.pi, 100)
    K=np.append(np.append(kxG_M,kxG_M), np.append(kxG_M,kxG_M))
    K=np.append(K, K)
    
    plt.plot(K, Edata, '^',label='0')
    plt.plot(K, engy(K,2.91,3.19,0.82,0.18), 'o',label='1')
    plt.show()
    popt, pcov = curve_fit(engy, K, Edata, maxfev=50000, p0=(2.91,3.19,0.82,0.18))
    
    #popt=[2.91,3.19,0.82, 0.18]
    perr = np.sqrt(np.diag(pcov))
    print(popt)
    print(perr)    
    return popt 
    
    
def findfit(i1):
    data=np.loadtxt("4/0-"+str('%.2f' % i1)+"-BAND.dat", delimiter=None)
    press = np.zeros(8, dtype=np.int)
    a=np.zeros((300, 8))
    index=0
    for i in range(0, (len(data)//300)):
        for j in range(100):
            if data[i*300+100+j, 1]*data[i*300+101+j, 1]<0:
                press[index]=i
                index=index+1
    for index in range(8):
        a[:, index]=data[press[index]*300:press[index]*300+300, 1]
    a[:, 1]=a[::-1, 1]    
    a[:, 3]=a[::-1, 3]    
    a[:, 5]=a[::-1, 5]   
    a[:, 7]=a[::-1, 7]    
    #x=np.linspace(1, 300, 300)
    #plt.plot(x, a[:, 0], '^',label='0')
    #plt.plot(x, a[:, 1], 'o',label='1')
    #plt.plot(x, a[:, 2], 'D',label='2')
    #plt.plot(x, a[:, 3], '*',label='3')
    #plt.plot(x, a[:, 4], '<',label='4')
    #plt.plot(x, a[:, 5], '>',label='5')
    #plt.legend(loc=4) 
    #plt.show()
    Edata = np.append(np.append(a[100:200, 0], a[100:200, 1]), np.append(a[100:200, 2], a[100:200, 3]))
    Edata = np.append(Edata, np.append(a[100:200, 4], a[100:200, 5]))
    Edata = np.append(Edata, np.append(a[100:200, 6], a[100:200, 7]))
    #print(Edata)
    tfit=fit(Edata)
    #tfit=[2.91,3.19,0.82, 0.18]
    return tfit
    
def mytry(i1):
    cha=i1
    data=np.loadtxt("BAND-1-4.dat", delimiter=None)
    tfit=[2.74,3.05,0.57,0.16]
    tfit=findfit(i1)
    hw=[]
    for iii in range(120):
        i=iii*2.5
        if i<100:
            kx = np.pi*(99-i)/99
            ky = 0
        elif i<200:
            kx = np.pi*(i-100)/99; ky = kx
        else:
            kx = np.pi; ky = np.pi*(299-i)/99
        hw.append(matr(tfit,kx,ky))
    x = np.linspace(1, 300, 300)
    x2 = np.linspace(1, 300, 120)
    real=np.reshape(data[:, 1], (300,-1), order='F')  
    for i1 in range(np.shape(real)[1]):
        if i1%2==1:
            real[:, i1]=real[::-1, i1]
    return x,real,x2,hw
    

def fig1(i1):
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(5,5))
    x,real,x2,hw=mytry(i1)
    plt.plot(x2, hw, '--')
    plt.plot(x,real)
    plt.axis([1,300,-4,4]) 
    plt.xticks([0, 100, 200, 300],['X','$\Gamma$','M','X'],fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.grid(axis="x")    
    plt.grid(axis="y")    
    plt.ylabel('$E(eV)$',fontsize='x-large')
    #plt.subplots_adjust(top=0.95, bottom=0.05, right=0.95, left=0.25, hspace=0.05, wspace=0.05)
    plt.margins(0, 0)
    plt.show()
    

def fig(i):
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(3,5))
    x,real,x2,hw=mytry(i)
    plt.plot(x2, hw, '--')
    plt.plot(x,real)
    plt.axis([1,300,-4,4]) 
    plt.xticks([0, 100, 200, 300],['X','$\Gamma$','M','X'],fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.grid(axis="x")    
    plt.grid(axis="y")    
    plt.margins(0, 0)
    plt.show()
    
fig(4)
