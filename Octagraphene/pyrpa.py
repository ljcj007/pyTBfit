import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import cmath

def matr(tpara,kx,ky,kz):
    t1 = tpara[0]; t2 = tpara[1]; t3 = tpara[2]; t4 = tpara[3]
    s=np.zeros((4,4), dtype=np.complex)
    s[0,1] = -1*t1; 
    s[0,2] = -1*t2*cmath.exp((0.+1.j)*ky)-t3;
    s[0,3] = -1*t1;
    s[1,2] = -1*t1;
    s[1,3] = -1*t2*cmath.exp((0.+1.j)*kx)-t3;
    s[2,3] = -1*t1;
    s[0,0] = -1*t4*cmath.cos(kz);
    s[1,1] = -1*t4*cmath.cos(kz);
    s[2,2] = -1*t4*cmath.cos(kz);
    s[3,3] = -1*t4*cmath.cos(kz);
    s=np.mat(s)
    return s+s.H
    
def eny4(tpara,kx,ky,kz):
    E,M=np.linalg.eig(matr(tpara,kx,ky,kz))
    return E.real,M
    
def diri(tpara,nk,nkz):
    eng=[]
    mat=[]
    for ix in range(nk):  
        kx=2*np.pi*ix/nk
        for iy in range(nk):
            ky=2*np.pi*iy/nk
            for iz in range(nkz):
                kz=2*np.pi*iz/nkz
                E,M=eny4(tpara,kx,ky,kz)
                eng.append(E)
                mat.append(M)
    return eng,mat

#def chi(tpara,nk,nkz,B):
#    eng,mat=diri(tpara,nk,nkz)
#    chi=[]
#    for q in range(1):  
#        qx=q//(nk*nkz);
#        qy=q%(nk*nkz)//nkz;
#        qz=q%(nk*nkz)%nkz;
#        chiq=np.zeros((16,16), dtype=np.complex)
#        for ik in range(nk*nk*nkz):  
#            kx=ik//(nk*nkz);
#            ky=ik%(nk*nkz)//nkz;
#            kz=ik%(nk*nkz)%nkz;
#            kx1=(kx+qx)%nk;
#            ky1=(ky+qy)%nk;
#            kz1=(kz+qz)%nkz;
#            ik1=kx1*nk*nkz+ky1*nkz+kz1
#            for iiii in range(16):
#                l1=iiii//4
#                l2=iiii%4
#                for jjjj in range(16):
#                    l3=jjjj//4
#                    l4=jjjj%4
#                    chi0=0
#                    for alpha in range(4):
#                        for beta in range(4):
#                            if (np.abs(eng[ik][alpha]-eng[ik1][beta])>1e-7):
#                                chi0=chi0+mat[ik1][l2,alpha]*np.conj(mat[ik1][l3,alpha])*mat[ik][l4,beta]\
#                                    *np.conj(mat[ik][l1,beta])\
#                                    *( 1/(np.exp(eng[ik1][beta]*B)+1)-1/(np.exp(eng[ik][alpha]*B)+1) )\
#                                    /(eng[ik][alpha]-eng[ik1][beta])
#                    chiq[iiii,jjjj]=chi0
#        chi.append(chiq)
#    return chi
def chi(tpara,nk,nkz,B):
    eng,mat=diri(tpara,nk,nkz)
    enp=( 1/(np.exp(eng*B)+1))
    chi=[]
    for ik1 in range(nk*nk*nkz):  
        chiq=np.zeros((16,16), dtype=np.complex)
        for iiii in range(16):
            l1=iiii//4
            l2=iiii%4
            for jjjj in range(16):
                l3=jjjj//4
                l4=jjjj%4
                chi0=0
                for ik in range(nk*nk*nkz):
                    for alpha in range(4):
                        for beta in range(4):
                            chi0=chi0+mat[ik1][l2,alpha]*np.conj(mat[ik1][l3,alpha])*mat[ik][l4,beta]\
                                 *np.conj(mat[ik][l1,beta])\
                                 *( 1/(np.exp(eng[ik1][beta]*B)+1)-1/(np.exp(eng[ik][alpha]*B)+1) )\
                                 /(eng[ik][alpha]-eng[ik1][beta]+0.00001j)
                chiq[iiii,jjjj]=chi0
        chi.append(chiq)
    return chi  
    
    
    

e,m=diri([2.7,2.9,0.4,0.2],1,1)
print(m[0])
e1=matr([2.7,2.9,0.4,0.2],0,0,0)
print(matr([2.7,2.9,0.4,0.2],0,0,0))
print(m[0].H*e1*m[0])
#chi=chi([2.7,2.9,0.4,0.2],2,2,1)
#print(chi[0])
#print(np.linalg.eig(np.mat(chi[0])))
