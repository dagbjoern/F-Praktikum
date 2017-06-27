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


def gerade(x,m,b):
    return m*x+b

def Helmholz_B(I,N,R):
    B=const.mu_0*(8*I*N/np.sqrt(125)*R)
    return B

def g_F(m):
    g=(4*np.pi*const.m_e)/(const.e*m)
    return g

Frequenz, drehung_1, off_1 , drehung_2 ,off_2 = np.genfromtxt('messwerte.txt',unpack=True)

Sweep_Spule_N=22
Sweep_Spule_R=0.1639

hori_Spule_N=154*2
hori_Spule_R=0.1579

vert_Spule_N=20*2
vert_Spule_R=0.11735

vertikal_und_swepp_ratio=0.1
hori_ratio=0.3

#strome Berechnen
I_1=drehung_1*vertikal_und_swepp_ratio
I_off_1=off_1*hori_ratio

I_2=drehung_2*vertikal_und_swepp_ratio
I_off_2=off_2*hori_ratio

#Magnetfelder berechnen
B_1_sweep=Helmholz_B(I_1,Sweep_Spule_N,Sweep_Spule_R)
B_1_horz=Helmholz_B(I_off_1,hori_Spule_N,hori_Spule_R)

B_2_sweep=Helmholz_B(I_2,Sweep_Spule_N,Sweep_Spule_R)
B_2_horz=Helmholz_B(I_off_2,hori_Spule_N,hori_Spule_R)


B_1_ges=B_1_sweep+B_1_horz
B_2_ges=B_2_sweep+B_2_horz

print(B_1_ges,B_2_ges)

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

b_1=unp.uarray(params_1[1],fehler_1[1])
b_2=unp.uarray(params_2[1],fehler_2[1])


g_F_1=g_F(m_1)
g_F_2=g_F(m_2)
print(g_F_1,g_F_2)
# print('B_1',b_1)
# print('B_2',b_2)



plt.figure(1)
plt.plot(Frequenz,B_1_ges*10**(6),'rx',label=r'Messwerte $Pik_1$')
plt.plot(Frequenz,B_2_ges*10**(6),'bx',label=r'Messwerte $Pik_2$')
plt.plot(Frequenz,gerade(Frequenz,*params_1)*10**(6),'-m',label=r'lin 1')
plt.plot(Frequenz,gerade(Frequenz,*params_2)*10**(6),'-c',label=r'lin 2')
plt.legend(loc='best')
plt.xlabel(r'Frequenz $\nu$ ')
plt.ylabel(r' Magnetischeflussdichte $B/ \mu T$')
plt.savefig('build/plot1_messwerte.pdf')
