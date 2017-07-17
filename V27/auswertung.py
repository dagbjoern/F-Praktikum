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

def lamDelta(lam,d,n):
    return lam**2/(2*d) *np.sqrt(1/(n**2-1))

def Auflosung(lam,L,n):
    return (L/lam) * (n**2-1)

def landefak(J,S,L):
    return(3*J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1))

def delta_E(m_1,g_1,m_2,g_2,B):
    mu_b=const.value("Bohr magneton")
    return g_12(m_1,g_1,m_2,g_2)*mu_b*B

def g_12(m_1,g_1,m_2,g_2):
    return (m_1*g_1-m_2*g_2)


I1 ,B1 =np.genfromtxt('messwertehyst1.txt',unpack=True)

I2 ,B2 =np.genfromtxt('messwertehyst2.txt',unpack=True)


lam_rot=643.8e-9
lam_blau=480e-9

#Maße der Lummer-Gehrcke-Platte
d=0.004
L=0.120
n_rot=1.4567
n_blau=1.4635


#lande-Faktoren
#1P1
J=1
S=0
L=1
g_1P1=landefak(J,S,L)
print('g_i von 1P1 =',g_1P1)

#1D2
S=0
L=2
J=2
g_1D2=landefak(J,S,L)
print('g_i von 1D2 =',g_1D2)

#3S1
S=1
L=0
J=1
g_3S1=landefak(J,S,L)
print('g_i von 3S1 =',g_3S1)


#3P1
S=1
L=1
J=1
g_3P1=landefak(J,S,L)
print('g_i von 3P1 =',g_3P1)

print('Despersionsgebiet rot=',lamDelta(lam_rot,d,n_rot))
print('Despersionsgebiet blau=',lamDelta(lam_blau,d,n_blau))
print('Auflöungsvermögen rot=',lamDelta(lam_rot,L,n_rot))
print('Auflöungsvermögen blau=',lamDelta(lam_blau,L,n_blau))

#rot
#sigma
#1P1 <-> 1D2
#def g_12(m_1,g_1,m_2,g_2):
print('rot')
print('landfaktor Delta m=-1', g_12(1,g_1P1,0,g_1D2))
print('landfaktor Delta m=1')


# blau
# 3S1  <-> 3P1
# sigma
print('blau')
print('landfaktor Delta m=-1', g_12(1,g_3S1,0,g_3P1))
print('landfaktor Delta m=-1', g_12(-1,g_3S1,0,g_3P1))





def hysterese(x,A,x0,C):
    return A*(np.arctan((x-x0)/C))
# pi
params_hyst , cov_hyst = curve_fit(hysterese ,I1,B1,p0=[1000,0,10])

fehler_hyst=np.sqrt(np.diag(cov_hyst))

A=unp.uarray(params_hyst[0],fehler_hyst[0])

x0=unp.uarray(params_hyst[1],fehler_hyst[1])

C=unp.uarray(params_hyst[2],fehler_hyst[2])

print('A=',A)
print('I_0=',x0)
print('C=',C)


#hysterese

I_lin=np.linspace(0,20)
plt.figure(1)
plt.plot(I1,B1,'xr',label=r'Messwerte')
plt.plot(I2,B2,'xr',)

#plt.plot(T1[55:],(I1[55:]),'bx',label=r'Untergrund')
plt.plot(I_lin,hysterese(I_lin,*params_hyst),'c-',label=r'Fit')
#plt.plot(I_lin,hysterese(I_lin,1000,0,10),'g-')
#plt.plot(x,gerade(x,m_untergrund1,b_untergrund1),'k-',label=r'Untergrund fit linregress')
# plt.plot(T1[:j],I1[:j],'mx',label=r'Messpunkte ohne Untergrund')
# plt.plot(T1[j:],I1[j:]-gerade(T1[j:],*params_untergund1),'mx',)
# plt.plot(T1,I1,'kx',label=r'Messpunkte mit Untergrund')
plt.legend(loc='best')
plt.xlabel(r'Strom $I/ \si{\ampere} $')
plt.ylabel(r'Magnetfeld $B/ \si{\milli\tesla}$')
plt.savefig('build/plot_hyst.pdf')
