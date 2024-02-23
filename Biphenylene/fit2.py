import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from cmath import exp
def eny4(tpara,kx,ky):
    #x1=0.193;y1=0.338;x2=0.500;y2=0.160
    x1=0.193;y1=0.337;x2=0.500;y2=0.158
    t1 = tpara[0]; t2 = tpara[1]; t3 = tpara[2]; t4 = tpara[3];
    t5 = tpara[4]; t6 = tpara[5]; t7 = tpara[6]; t8 = tpara[7]
    t9 = tpara[8]; t10 = tpara[9]; t11 = tpara[10]; t12 = tpara[11]
    t13 = tpara[12]; t14 = tpara[13]
    #t1 = tpara[0];t2=t1*1.455/1.402;t3=t1*1.455/1.439; t4 = t1*1.455/1.449
    s=np.zeros((6,6), dtype=complex)
    s[0,1] = -1*(t1)*exp((0.+1.j)*(ky*(1-2*y1)))-t12*exp((0.+1.j)*(ky*(-2*y1))); #t1=1.455
    s[1,2] = -1*(t2)*exp((0.+1.j)*(kx*(x2-x1)+ky*(y1-y2)))-t10*exp((0.+1.j)*(kx*(-x2-x1)+ky*(y1-y2)));#t2=1.402
    s[2,3] = -1*(t2)*exp((0.+1.j)*(kx*(x2-x1)+ky*(y2-y1)))-t10*exp((0.+1.j)*(kx*(-x2-x1)+ky*(y2-y1)));
    s[3,4] = -1*(t1)*exp((0.+1.j)*(ky*(2*y1-1)))-t12*exp((0.+1.j)*(ky*(2*y1))); #t12=3.049
    s[4,5] = -1*(t2)*exp((0.+1.j)*(kx*(x1-x2)+ky*(y2-y1)))-t10*exp((0.+1.j)*(kx*(x1+x2)+ky*(y2-y1)));
    s[5,0] = -1*(t2)*exp((0.+1.j)*(kx*(x1-x2)+ky*(y1-y2)))-t10*exp((0.+1.j)*(kx*(x1+x2)+ky*(y1-y2)));#t10=2.720
    s[0,4] = -1*(t4)*exp((0.+1.j)*(kx*(-2*x1)))-t6*exp((0.+1.j)*(kx*(2*x2-2*x1))); #t4=1.449,
    s[1,3] = -1*(t4)*exp((0.+1.j)*(kx*(-2*x1)))-t6*exp((0.+1.j)*(kx*(2*x2-2*x1))); #t6=2.298,
    s[2,5] = -1*(t3)*exp((0.+1.j)*(ky*(2*y2)))-(t11)*exp((0.+1.j)*(ky*(2*y2-1))); #t3=1.439,t11=3.061
    s[0,3] = -1*(t5)*exp((0.+1.j)*(kx*(-2*x1)+ky*(1-2*y1)))-t9*exp((0.+1.j)*(kx*(2*x2-2*x1)+ky*(1-2*y1)))-t14*exp((0.+1.j)*(kx*(-2*x1)+ky*(-2*y1))); #t5=2.053,t14=3.376
    s[1,4] = -1*(t5)*exp((0.+1.j)*(kx*(-2*x1)+ky*(2*y1-1)))-t9*exp((0.+1.j)*(kx*(2*x1-2*x2)+ky*(2*y1-1)))-t14*exp((0.+1.j)*(kx*(-2*x1)+ky*(-2*y1))); #t9=2.720
    s[0,2] = -1*t7*exp((0.+1.j)*(kx*(x2-x1)+ky*(1-y2-y1)))-t8*exp((0.+1.j)*(kx*(x2-x1)+ky*(-y2-y1)))-t13*exp((0.+1.j)*(kx*(-x1-x2)+ky*(-y1-y2)));#t7=2.534
    s[2,4] = -1*t7*exp((0.+1.j)*(kx*(x2-x1)+ky*(y2+y1-1)))-t8*exp((0.+1.j)*(kx*(x2-x1)+ky*(y2+y1)))-t13*exp((0.+1.j)*(kx*(-x1-x2)+ky*(y1+y2)));#t13=3.434
    s[1,5] = -1*t7*exp((0.+1.j)*(kx*(x2-x1)+ky*(y2+y1-1)))-t8*exp((0.+1.j)*(kx*(x2-x1)+ky*(y2+y1)))-t13*exp((0.+1.j)*(kx*(-x1-x2)+ky*(y1+y2)));#t8=2.521
    s[3,5] = -1*t7*exp((0.+1.j)*(kx*(x1-x2)+ky*(y2+y1-1)))-t8*exp((0.+1.j)*(kx*(x1-x2)+ky*(y2+y1)))-t13*exp((0.+1.j)*(kx*(x1+x2)+ky*(y1+y2)));
    s=np.mat(s)
    s=s+s.H
    E,_=np.linalg.eig(s)
    t=E.real
    return sorted(t)
      
def eny1(tpara,kx,ky):
    x1=0.193
    y1=0.338
    x2=0.500
    y2=0.160
    t1 = tpara[0]; t2 = tpara[1]; t3 = tpara[2]; t4 = tpara[3];
    t5 = tpara[4]; t6 = tpara[5]; t7 = tpara[6]; t8 = tpara[7]
    #t9 = tpara[4]; t10 = tpara[5]; t11 = tpara[6]; t12 = tpara[7]
    #t1 = tpara[0];t2=t1*1.455/1.402;t3=t1*1.455/1.439; t4 = tpara[3]
    s=np.zeros((6,6), dtype=complex)
    s[0,1] = -1*(t1)*exp((0.+1.j)*(ky*(1-2*y1)))#-t8*exp((0.+1.j)*(ky*(-2*y2-2*y1))); #t1=1.455
    s[1,2] = -1*(t2)*exp((0.+1.j)*(kx*(x2-x1)+ky*(y1-y2)))-t7*exp((0.+1.j)*(kx*(-x2-x1)+ky*(y1-y2)));#t2=1.402
    s[2,3] = -1*(t2)*exp((0.+1.j)*(kx*(x2-x1)+ky*(y2-y1)))-t7*exp((0.+1.j)*(kx*(-x2-x1)+ky*(y2-y1)));
    s[3,4] = -1*(t1)*exp((0.+1.j)*(ky*(2*y1-1)))#-t8*exp((0.+1.j)*(ky*(2*y2+2*y1))); #t8=3.049
    s[4,5] = -1*(t2)*exp((0.+1.j)*(kx*(x1-x2)+ky*(y2-y1)))-t7*exp((0.+1.j)*(kx*(x1+x2)+ky*(y2-y1)));
    s[5,0] = -1*(t2)*exp((0.+1.j)*(kx*(x1-x2)+ky*(y1-y2)))-t7*exp((0.+1.j)*(kx*(x1+x2)+ky*(y1-y2)));#t7=2.720
    s[0,4] = -1*(t1)*exp((0.+1.j)*(kx*(-2*x1)))-t5*exp((0.+1.j)*(kx*(2*x2-2*x1))); #t1=1.449,
    s[1,3] = -1*(t1)*exp((0.+1.j)*(kx*(-2*x1)))-t5*exp((0.+1.j)*(kx*(2*x2-2*x1))); #t5=2.298,
    s[2,5] = -1*(t3)*exp((0.+1.j)*(ky*(2*y2)))#-(t8)*exp((0.+1.j)*(ky*(2*y2-1))); #t3=1.439,t8=3.061
    s[0,3] = -1*(t4)*exp((0.+1.j)*(kx*(-2*x1)+ky*(1-2*y1)))#-t7*exp((0.+1.j)*(kx*(2*x2-2*x1)+ky*(1-2*y1))); #t4=2.053
    s[1,4] = -1*(t4)*exp((0.+1.j)*(kx*(-2*x1)+ky*(2*y1-1)))#-t7*exp((0.+1.j)*(kx*(2*x1-2*x2)+ky*(2*y1-1))); #t7=2.720
    s[0,2] = -1*t6*exp((0.+1.j)*(kx*(x2-x1)+ky*(1-y2-y1)))-t8*exp((0.+1.j)*(kx*(x2-x1)+ky*(y2+y1-1)));#t6=2.534
    s[2,4] = -1*t6*exp((0.+1.j)*(kx*(x2-x1)+ky*(y2+y1-1)))-t8*exp((0.+1.j)*(kx*(x2-x1)+ky*(1-y2-y1)));
    s[1,5] = -1*t6*exp((0.+1.j)*(kx*(x2-x1)+ky*(y2+y1-1)))-t8*exp((0.+1.j)*(kx*(x2-x1)+ky*(1-y2-y1)));#t6=2.521
    s[3,5] = -1*t6*exp((0.+1.j)*(kx*(x1-x2)+ky*(y2+y1-1)))-t8*exp((0.+1.j)*(kx*(x2-x1)+ky*(1-y2-y1)));
    s=np.mat(s)
    s=s+s.H
    E,_=np.linalg.eig(s)
    t=E.real
    return sorted(t)
    
    
def engy(K, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14):
    nk=len(K)
    an=np.zeros(nk)
    tpara=[t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14]
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
    #Edata = np.append(real[:, 15], real[:, 16])
    kxG_M = np.linspace(0,500,500)
    popt, pcov = curve_fit(engy, kxG_M, Edata, maxfev=50000, p0=(2.35143193,  3.31765141,  4.42764936,  2.92444521, -0.31549675,  0.66039258,
 -0.34997734,  0.61737266, -0.26438026, -0.52132035,  0.82047972,  0.23464918,  0.4245566,   0.22497889))
    #popt, pcov = curve_fit(engy, kxG_M, Edata, maxfev=50000, p0=(2.9,3.1,2.9,2.9,1.2,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.0,0.0))
    perr = np.sqrt(np.diag(pcov))
    print(popt)
    print(perr)
    return popt,perr

def engy1(K, t1, t2, t3, t4, t5, t6, t7, t8):
    nk=len(K)
    an=np.zeros(nk)
    tpara=[t1, t2, t3, t4, t5, t6, t7, t8]
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
    popt, pcov = curve_fit(engy1, kxG_M, Edata, maxfev=50000, p0=(2.7,2.9,2.7,2.7,0.0,0.0,0.0,0.0))
    perr = np.sqrt(np.diag(pcov))
    print(popt)
    print(perr)
    return popt,perr


def engy2(K, t1, t2, t3, t4, t5, t6, t7, t8):
    nk=len(K)
    an=np.zeros(nk)
    tpara=[t1, t2, t3, t4, t5, t6, t7, t8]
    for i in range(0, nk-1):
        if i<50:
            kx = np.pi; ky = np.pi*(i)/49
            t=eny4(tpara,kx,ky)
            an[i]=t[2]
        elif i<100:
            kx = np.pi; ky = np.pi*(i-50)/49
            t=eny4(tpara,kx,ky)
            an[i]=t[3]
    return an

def findfit2(data):    
    real=np.reshape(data[:, 1], (250,-1), order='F')  
    for i1 in range(np.shape(real)[1]):
        if i1%2==1:
            real[:, i1]=real[::-1, i1]
    Edata = np.append(real[50:100, 11], real[50:100, 12])
    kxG_M = np.linspace(0,100,100)
    popt, pcov = curve_fit(engy2, kxG_M, Edata, maxfev=50000, p0=(2.7,2.7,2.7,2.7,0.0,0.0,0.0,0.0))
    perr = np.sqrt(np.diag(pcov))
    print(popt)
    print(perr)
    return popt,perr

def energy_bond(data):
    tfit,perr=findfit(data)
    #tfit=[2.71836897, 3.01048793, 2.71604796, 0.70721011]
    #tfit=[2.94292,  2.92715985,  -0.56,1.29]
    #tfit=[2.87443588, 3.06050229, 2.66715071, 3.03441115]
    #tfit=[3.0,  2.4,  2.3, 3.6]
    #tfit=[2.79489393,  3.28674952, 3.13478288,  2.88922255]
     #[ 2.26513004  3.18340131  2.94907984  2.521327    0.08634244  0.14394721
    #tfit=[2.2513004,  3.18340131,  2.94907984,  2.521327,    0.08634244,  0.14394721, -0.21099235,  0.21867582,  0.18220586, -0.2747581,   1.09864612,  0.13913331] 
    #tfit=[2.2513004,  3.18340131,  2.94907984,  2.521327,    0.08634244,  0.14394721, -0.21099235,  0.21867582,  0.18220586, -0.2747581,   1.09864612,  0.13913331,0,0]     
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
    data=np.loadtxt("./BAND-%.2d.dat" % i, delimiter=None)
    x = np.linspace(1, 250, 250)
    real=np.reshape(data[:, 1], (250,-1), order='F')  
    for i1 in range(np.shape(real)[1]):
        if i1%2==1:
            real[:, i1]=real[::-1, i1]
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(3,4))
    plt.plot(x2, hw, linestyle='--')#-3.65-0.59
    plt.plot(x,real)
    plt.axis([1,250,-5,5]) 
    plt.xticks([0, 50, 100, 150, 200,250],['$\Gamma$','K','M','$\Gamma$','K`','M'],fontsize='x-large')
    plt.yticks([-4, -2, 0, 2, 4],fontsize='x-large')
    plt.grid(axis="x")    
    plt.grid(axis="y")    
    plt.margins(0, 0)
    plt.subplots_adjust(top=0.98, bottom=0.10, right=0.98, left=0.12, hspace=0.05, wspace=0.05)
    plt.savefig('./'+str('%.2d' % i)+'.jpg',dpi=600,format='jpg')
    plt.show()
fig(11)

