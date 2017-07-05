import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const
from uncertainties import ufloat
from scipy.integrate import trapz
from scipy.integrate import simps
from tabulate import tabulate

def abweichung(theo,mess):
    return(  (theo-mess)/theo)

def gerade(x,m,b):
    return m*x+b

def Helmholz_B(I,N,R):
    B=const.mu_0*(8*I*N/(np.sqrt(125)*R))
    return B

def g_F(m):
    g=(4*np.pi*const.m_e)/(const.e*m)
    return g

def I_spin(g_f,J):
    g_j=2.0023
    return J*((g_j/g_f) -1)

Frequenz, drehung_1, off_1 , drehung_2 ,off_2 = np.genfromtxt('messwerte.txt',unpack=True)
Frequenz=Frequenz*1000
drehung_vertikal=2.2

Sweep_Spule_N=11
Sweep_Spule_R=0.1639

hori_Spule_N=154
hori_Spule_R=0.1579

vert_Spule_N=20
vert_Spule_R=0.11735

vertikal_und_swepp_ratio=0.1
hori_ratio=0.3


I_vertikal=drehung_vertikal*vertikal_und_swepp_ratio
B_vertikal=Helmholz_B(I_vertikal,vert_Spule_N,vert_Spule_R)

print('Strom vertikal=', I_vertikal,'\nFlussdichte B_vertikal=',B_vertikal )
#strome Berechnen
I_1=drehung_1*vertikal_und_swepp_ratio
I_off_1=off_1*hori_ratio
print('\nI_1',I_1,'\n Ioff_1',I_off_1)

I_2=drehung_2*vertikal_und_swepp_ratio
I_off_2=off_2*hori_ratio
print('\nI_2',I_2,'\nIoff_2',I_off_2)


#Magnetfelder berechnen
B_1_sweep=Helmholz_B(I_1,Sweep_Spule_N,Sweep_Spule_R)
B_1_horz=Helmholz_B(I_off_1,hori_Spule_N,hori_Spule_R)

B_2_sweep=Helmholz_B(I_2,Sweep_Spule_N,Sweep_Spule_R)
B_2_horz=Helmholz_B(I_off_2,hori_Spule_N,hori_Spule_R)


B_1_ges=B_1_sweep+B_1_horz
B_2_ges=B_2_sweep+B_2_horz

print('B_1ges',B_1_ges,'\nB_2ges',B_2_ges)

#
# def fehler_b(std_m,x):
#     return(std_m*np.sqrt(np.mean(x**2)))


params_1 , cov_1 = curve_fit(gerade ,Frequenz,B_1_ges)
params_2 , cov_2 = curve_fit(gerade ,Frequenz,B_2_ges)

# print('lin ',std_m_untergrund1,std_b_untergrund1)
# print('params',params_untergund1)
# print('fehler',np.sqrt(np.diag(cov_untergrund1)))
fehler_1=np.sqrt(np.diag(cov_1))
fehler_2=np.sqrt(np.diag(cov_2))

m_1=unp.uarray(params_1[0],fehler_1[0])
m_2=unp.uarray(params_2[0],fehler_2[0])

print('m_1',m_1)
print('m_2',m_2)
b_1=unp.uarray(params_1[1],fehler_1[1])
b_2=unp.uarray(params_2[1],fehler_2[1])
Erde_hor=unp.uarray([params_1[1],params_2[1]],[fehler_1[1],fehler_2[1]])
print()
print('b_1',b_1)
print('b_2',b_2)
print(Erde_hor.mean())
#print(unp.mean([b_1,b_2]))


g_F_1=g_F(m_1)
g_F_2=g_F(m_2)
print('g_1',g_F_1)
print('g_2',g_F_2)
# print('B_1',b_1)
# print('B_2',b_2)
J=0.5
print('I_1',I_spin(g_F_1,J))
print('I_2',I_spin(g_F_2,J))

x=np.linspace(-1,np.amax(Frequenz)+100000)
plt.figure(1)
plt.plot(Frequenz/1000000,B_1_ges*10**(6),'rx',label=r'Messwerte $Pik_1$')
plt.plot(Frequenz/1000000,B_2_ges*10**(6),'bx',label=r'Messwerte $Pik_2$')
plt.plot(x/1000000,gerade(x,*params_1)*10**(6),'-r',alpha=0.25,label=r'Fit 1')
plt.plot(x/1000000,gerade(x,*params_2)*10**(6),'-b',alpha=0.25,label=r'Fit 2')
plt.legend(loc='best')
plt.xlabel(r'Frequenz $\nu$ / MHz ')
plt.ylabel(r' Magnetischeflussdichte $B/ \mu T$')
plt.grid('on')
plt.savefig('build/plot1_messwerte.pdf')

print('energie',1e6*const.h)


#quadratischer Zeemaneffekt
def zee(g_f,B,M_F,dE) :
    print('mu^2',(const.mu_0*const.mu_0))
    print('g_f',(g_f*g_f))
    print('B^2',B*B)
    return (g_f*g_f)*(const.mu_0*const.mu_0)*(B*B)*((1-2*M_F)/dE)


dE_87=4.54*10**(-24)
dE_85=2.01*10**(-24)
print(dE_85)
g_87=g_F_1
g_85=g_F_2
print(np.amax(B_1_ges))
print(np.amax(B_2_ges))
print('zee^2 E_87=',zee(g_87,np.amax(B_1_ges),1,dE_87) )
print('zee^2 E_85=',zee(g_85,np.amax(B_2_ges),1,dE_85) )


print('erde vertikal',abweichung(45.2e-6,B_vertikal))
print(19.3+0.891)
print('erde horiz',abweichung((19.3+0.891)*1e-6,Erde_hor.mean()))
print('abweichung I_Rb_87', abweichung(3/2,I_spin(g_F_1,J)))
print('abweichung I_Rb_85', abweichung(5/2,I_spin(g_F_2,J)))
print('abweichung verhältnis', abweichung(0.39,0.65))
