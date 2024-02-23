import matplotlib.pyplot as plt
import numpy as np



def energy_bond(data):
    x = np.linspace(1, 250, 250)
    #x=np.append(x,y)
    real=np.reshape(data[:, 1], (250,-1), order='F')  
    for i1 in range(np.shape(real)[1]):
        if i1%2==1:
            real[:, i1]=real[::-1, i1]
    return x,real

def figxx(i):
    data=np.loadtxt("./BAND-%.2d.dat" % i, delimiter=None)
    x,real=energy_bond(data)
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(5,5))
    plt.plot(x,real)
    plt.axis([1,250,-35,15]) 
    plt.xticks([0, 50, 100, 150, 200, 250],['$\Gamma$','K','M','$\Gamma$','K`','M'],fontsize='x-large')
    #plt.yticks([-1, -0.5, 0, 0.5, 1],fontsize='x-large')
    plt.grid(axis="x")    
    plt.grid(axis="y")    
    plt.margins(0, 0)
    plt.subplots_adjust(top=0.98, bottom=0.10, right=0.98, left=0.12, hspace=0.05, wspace=0.05)
    plt.savefig('./'+str('%.2d' % i)+'.eps',dpi=600,format='eps')
    plt.show()
    
    
def figx(i):
    data=np.loadtxt("./BAND-%.2d.dat" % i, delimiter=None)
    x,real=energy_bond(data)
    plt.rc('font',family='Times New Roman')
    plt.figure(figsize=(5,5))
    plt.plot(x[:],real[:,:])
    plt.axis([1,250,-5,5]) 
    plt.xticks([0, 50, 100, 150, 200, 250],['$\Gamma$','K','M','$\Gamma$','K`','M'],fontsize='x-large')
    #plt.yticks([-1, -0.5, 0, 0.5, 1],fontsize='x-large')
    plt.grid(axis="x")    
    plt.grid(axis="y")    
    plt.margins(0, 0)
    plt.subplots_adjust(top=0.98, bottom=0.10, right=0.98, left=0.12, hspace=0.05, wspace=0.05)
    plt.savefig('./'+str('%.2d' % i)+'.eps',dpi=600,format='eps')
    plt.show()
figx(11)
