import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import cmath
import math

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
    E,_=np.linalg.eig(np.hstack((np.vstack((s, s1)), np.vstack((s1, s)))))
    t=E.real
    return sorted(t)
    
def matr3(tpara,kx,ky):
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
    v1=np.vstack((np.vstack((s, s1)), s0))
    v2=np.vstack((np.vstack((s1, s)), s1))
    v3=np.vstack((np.vstack((s0, s1)), s))
    E,_=np.linalg.eig(np.hstack((np.hstack((v1, v2)), v3))) 
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
            an[i] = t[2]
        elif i<200:
            an[i] = t[3]
        elif i<300:
            an[i] = t[4]
        elif i<400:
            an[i] = t[5]
        else:
            an[i] = t[0]
    return an 
    
def fit(Edata):
    kxG_M = np.linspace(0, np.pi, 100)
    K=np.append(np.append(kxG_M,kxG_M), np.append(kxG_M,kxG_M))
    popt, pcov = curve_fit(engy, K, Edata, maxfev=50000, p0=(2.91,3.19,0.82, 0.23))
    perr = np.sqrt(np.diag(pcov))
    print(popt)
    print(perr)    
    return popt 
    
    
def findfit(i1):
    data=np.loadtxt("2/0-"+str('%.2f' % i1)+"-BAND.dat", delimiter=None)
    press = np.zeros(4, dtype=np.int)
    a=0
    for i in range(0, (len(data)//300)):
        for j in range(100):
            if data[i*300+100+j, 1]*data[i*300+101+j, 1]<0:
                press[a]=i
                a=a+1
    #x=np.linspace(1, 300, 300)
    a0=data[press[0]*300:press[0]*300+300, 1]
    a1=data[press[1]*300:press[1]*300+300, 1]
    a1=a1[::-1]    
    a2=data[press[2]*300:press[2]*300+300, 1]
    a3=data[press[3]*300:press[3]*300+300, 1]
    a3=a3[::-1]
    #plt.plot(x, a0, '^',label='0')
    #plt.plot(x, a1, '*',label='1')
    #plt.plot(x, a2, 'D',label='2')
    #plt.plot(x, a3, 'o',label='3')
    #plt.show()
    Edata = np.append(np.append(a0[100:200], a1[100:200]), np.append(a2[100:200], a3[100:200]))
    #tfit=[-1.13224772e+01, -2.13484988e+00 , 2.36586008e-01, 5.81806723e-09]
    #print(sorted(tfit))
    #print(Edata)
    tfit=fit(Edata)
    return tfit
    
def mytry(i1):
    cha=i1
    data=np.loadtxt("BAND-1-2.dat", delimiter=None)
    tfit=[2.74,3.05,0.57,0.16]
    tfit=findfit(2)
    #tfit=[2.70,2.96,0.52,0.18]
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
        hw.append(matr(tfit,kx,ky))
    x = np.linspace(1, 300, 300)
    x2 = np.linspace(1, 300, 60)
    real=np.reshape(data[:, 1], (300,-1), order='F')  
    for i1 in range(np.shape(real)[1]):
        if i1%2==1:
            real[:, i1]=real[::-1, i1]
    return x,real,x2,hw

def findFS():
    tfit=findfit(2)
    hwx=[]
    hwy=[]
    hqx=[]
    hqy=[]
    eng1=0
    eng2=0
    en1=0
    en2=0
    eng3=0
    eng4=0
    en3=0
    en4=0
    for i in range(201):
        for j in range(201):
            kx = np.pi*(i-100)/100
            ky = np.pi*(j-100)/100
            en1=matr(tfit,kx,ky)[2]
            en2=matr(tfit,kx,ky)[3]
            en3=matr(tfit,kx,ky)[4]
            en4=matr(tfit,kx,ky)[5]
            if en1*eng1<0:
                hwx.append(kx)
                hwy.append(ky)
            if en3*eng3<0:
                hqx.append(kx)
                hqy.append(ky)
            if en2*eng2<0:
                hwx.append(kx)
                hwy.append(ky)
            if en4*eng4<0:
                hqx.append(kx)
                hqy.append(ky)
            eng1=en1
            eng2=en2
            eng3=en3
            eng4=en4
    eng1=0
    eng2=0
    eng3=0
    eng4=0
    for j in range(201):
        for i in range(201):
            kx = np.pi*(i-100)/100
            ky = np.pi*(j-100)/100
            en1=matr(tfit,kx,ky)[2]
            en2=matr(tfit,kx,ky)[3]
            en3=matr(tfit,kx,ky)[4]
            en4=matr(tfit,kx,ky)[5]
            if en1*eng1<0:
                hwx.append(kx)
                hwy.append(ky)
            if en3*eng3<0:
                hqx.append(kx)
                hqy.append(ky)
            if en2*eng2<0:
                hwx.append(kx)
                hwy.append(ky)
            if en4*eng4<0:
                hqx.append(kx)
                hqy.append(ky)
            eng1=en1
            eng2=en2
            eng3=en3
            eng4=en4
    return hwx, hwy, hqx, hqy

def figFM():
    plt.figure(figsize=(3,3))
    hwx, hwy, hqx, hqy=findFS()
    x= np.linspace(0, math.pi, 100)
    y1=0*x
    y2=x
    y3=math.pi+y1
    plt.xlim(-math.pi, math.pi)
    plt.ylim(-math.pi, math.pi)
    plt.plot(x, y1, 'b')
    plt.plot(x, y2, 'b')
    plt.plot(y3, x, 'b')
    plt.scatter(y1,y1,s=2)
    plt.scatter(hwx,hwy,s=2)
    plt.scatter(hqx,hqy,s=2)
    plt.xticks([-math.pi, 0, math.pi],['$-\pi$','0','$\pi$'], fontsize='x-large')
    plt.yticks([-math.pi, 0, math.pi],['$-\pi$','0','$\pi$'], fontsize='x-large')
    plt.xlabel('$k_x$', fontsize='x-large')
    plt.xlabel('$k_y$', fontsize='x-large')
    plt.annotate("$\Gamma$" , xy=(0,0), xytext=(1, 0), textcoords='offset points', fontsize='x-large')
    plt.annotate("$M$" , xy=(math.pi,math.pi), xytext=(1, 0), textcoords='offset points', fontsize='x-large')
    plt.annotate("$X$" , xy=(math.pi, 0), xytext=(1, 0), textcoords='offset points', fontsize='x-large')
    plt.subplots_adjust(top=0.95, bottom=0.20, right=0.95, left=0.20)
    plt.savefig('./femi2.eps',dpi=600,format='eps')
    plt.show()

def fig1(i1):
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(3,4))
    x,real,x2,hw=mytry(i1)
    plt.plot(x,real)
    plt.axis([1,300,-4.5,4.5]) 
    plt.xticks([0, 100, 200, 300],['X','$\Gamma$','M','X'],fontsize='x-large')
    plt.yticks([-4, -2, 0, 2, 4],fontsize='x-large')
    plt.grid(axis="x")    
    plt.grid(axis="y")    
    plt.subplots_adjust(top=0.98, bottom=0.10, right=0.98, left=0.12, hspace=0.05, wspace=0.05)
    plt.savefig('./'+str('%.2d' % i1)+'.eps',dpi=600,format='eps')
    plt.show()
    
def fig():
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(3,3))
    x,real,x2,hw=mytry(2)
    plt.plot(x2, hw, 'o', alpha=0.5)
    plt.plot(x,real)
    plt.axis([170,230,-4.2,-0.8]) 
    plt.xticks([200], ['M'], fontsize='x-large')
    plt.yticks([-4, -3,-2, -1],fontsize='x-large')
    plt.grid(axis="x")      
    plt.ylabel('Energy (eV)',fontsize='x-large')
    plt.subplots_adjust(top=0.98, bottom=0.10, right=0.98, left=0.12, hspace=0.05, wspace=0.05)
    plt.margins(0, 0)
    plt.savefig('Mpoint2.eps',dpi=600,format='eps')
    plt.show()
    
def figt():
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(3,3))
    x,real,x2,hw=mytry(2)
    i=2
    plt.plot(x2, hw, linestyle='--')
    plt.plot(x,real[:,:])
    plt.axis([1,300,-4.5,4.5]) 
    plt.xticks([0, 100, 200, 300],['X','$\Gamma$','M','X'],fontsize='x-large')
    plt.yticks([-4, -2, 0, 2, 4],fontsize='x-large')
    plt.grid(axis="x")    
    plt.grid(axis="y")    
    plt.margins(0, 0)
    plt.text(130, 3.5, r'$\epsilon_1$',fontsize = 'xx-large') 
    plt.text(170, 0.5, r'$\epsilon_2$',fontsize = 'xx-large') 
    plt.text(110, -0.8, r'$\epsilon_3$',fontsize = 'xx-large') 
    plt.text(150, -3.4, r'$\epsilon_4$',fontsize = 'xx-large')
    plt.ylabel('Energy (eV)',fontsize='x-large') 
    plt.subplots_adjust(top=0.98, bottom=0.10, right=0.98, left=0.12, hspace=0.05, wspace=0.05)
    plt.savefig('./'+str('%.2d' % i)+'.eps',dpi=600,format='eps')
    plt.show()
    
figt()
