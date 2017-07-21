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

DlamD_rot=lamDelta(lam_rot,d,n_rot)
DlamD_blau=lamDelta(lam_blau,d,n_blau)
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
plt.plot(I_lin,hysterese(I_lin,*params_hyst),'c-',label=r'Fit-Funktion')
#plt.plot(I_lin,hysterese(I_lin,1000,0,10),'g-')
#plt.plot(x,gerade(x,m_untergrund1,b_untergrund1),'k-',label=r'Untergrund fit linregress')
# plt.plot(T1[:j],I1[:j],'mx',label=r'Messpunkte ohne Untergrund')
# plt.plot(T1[j:],I1[j:]-gerade(T1[j:],*params_untergund1),'mx',)
# plt.plot(T1,I1,'kx',label=r'Messpunkte mit Untergrund')
plt.legend(loc='best')
plt.xlabel(r'Strom $I/ \si{\ampere} $')
plt.ylabel(r'Magnetfeld $B/ \si{\milli\tesla}$')
plt.savefig('build/plot_hyst.pdf')


rot_Ds= 1510-1280  # D entspricht Delta
print(rot_Ds)
rot_ds1=1360-1220  # d entspricht delta
print(rot_ds1)
rot_ds2=1555-1440
print(rot_ds2)

rot_ds=unp.uarray(np.mean([rot_ds1,rot_ds2]),np.std([rot_ds1,rot_ds2]))
print(rot_ds)
blau_sig_Ds=2725-2515
print(blau_sig_Ds)
blau_sig_ds1=2580-2480
print(blau_sig_ds1)
blau_sig_ds2=2780-2685
print(blau_sig_ds2)

blau_sig_ds=unp.uarray(np.mean([blau_sig_ds1,blau_sig_ds2]),np.std([blau_sig_ds1,blau_sig_ds2]))
print(blau_sig_ds)

blau_pi_Ds=2945-2740
print(blau_pi_Ds)
blau_pi_ds1=2790-2710
print(blau_pi_ds1)
blau_pi_ds2=2985-2905
print(blau_pi_ds2)
blau_pi_ds=unp.uarray(np.mean([blau_pi_ds1,blau_pi_ds2]),np.std([blau_pi_ds1,blau_pi_ds2]))
print(blau_pi_ds)


def d_lam(Ds,ds,DlamD):
    return (1/2) * (ds/Ds)*DlamD

def d_E_lam(lam,d_lam):
    return -const.h*(const.c/lam**2)*d_lam

def g_mess(d_E,B):
    mu_b=const.value("Bohr magneton")
    return d_E/( mu_b*B)

def hysterese_unp(x,A,x0,C):
    return A*(unp.arctan((x-x0)/C))

print('A=',A)
print('I_0=',x0)
print('C=',C)

I_rot=unp.uarray(11.8,0.1)
I_blau_sig=unp.uarray(5.1,0.1)
I_blau_pi=unp.uarray(17.0,0.1)


B_rot=hysterese_unp(I_rot,A,x0,C)/1000
B_blau_sig=hysterese_unp(I_blau_sig,A,x0,C)/1000
B_blau_pi=hysterese_unp(I_blau_pi,A,x0,C)/1000

print('rot B',B_rot)
d_lam_rot=d_lam(rot_Ds,rot_ds,DlamD_rot)
print('d_lam-rot',d_lam_rot)
d_E_lam_rot=d_E_lam(lam_rot,d_lam_rot)
print('d_E_lam_rot',d_E_lam_rot/const.e)
print('g_rot',g_mess(d_E_lam_rot,B_rot))

g_ij_rot=g_mess(d_E_lam_rot,B_rot)

print('blau B_sig',B_blau_sig)
d_lam_sig_blau=d_lam(blau_sig_Ds,blau_sig_ds,DlamD_blau)
print('d_lam_sig_blau',d_lam_sig_blau)
d_E_lam_blau_sig=d_E_lam(lam_blau,d_lam_sig_blau)
print('d_E_lam_blau_sig',d_E_lam_blau_sig/const.e)
print('g_blau_sig',g_mess(d_E_lam_blau_sig,B_blau_sig))

g_ij_blau_sig=g_mess(d_E_lam_blau_sig,B_blau_sig)

print('blau_B_pi',B_blau_pi)
d_lam_pi_blau=d_lam(blau_pi_Ds,blau_pi_ds,DlamD_blau)
print('d_lam_pi_blau',d_lam_pi_blau)
d_E_lam_blau_pi=d_E_lam(lam_blau,d_lam_pi_blau)
print('d_E_lam_blau_pi',d_E_lam_blau_pi/const.e)
print('g_blau_pi',g_mess(d_E_lam_blau_pi,B_blau_pi))

g_ij_blau_pi=g_mess(d_E_lam_blau_pi,B_blau_pi)

def abweich(mess,theo):
    return(mess-theo)/theo

#Abweichungen
print('abweich rot sigma',abweich(np.abs(g_ij_rot),1),'theorie',1)
print('abweich blau sigma',abweich(np.abs(g_ij_blau_sig),(2+1.5)/2),'theorie',(2+1.5)/2)
print('abweich blau pi',abweich(np.abs(g_ij_blau_pi),0.5),'theorie',0.5)
