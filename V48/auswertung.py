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

I1, T1 = np.genfromtxt('messung1.txt',unpack=True)

I2, T2 = np.genfromtxt('messung2.txt',unpack=True)

T1=T1+273.15
T2=T2+273.15


#
# def fehler_b(std_m,x):
#     return(std_m*np.sqrt(np.mean(x**2)))


def gerade(x,m,b):
    return m*x+b

def ef(x,a,b,c):
    return a*np.exp(c*x)+b
#linerare regression

# m_untergrund1 , b_untergrund1 , r ,p ,std_m_untergrund1 =stats.linregress(T1[55:],(I1[55:]))

# std_b_untergrund1=fehler_b(std_m_untergrund1,T1[55:])

params_untergund1 , cov_untergrund1 = curve_fit(gerade ,(T1[55:]),(I1[55:]))


# print('lin ',std_m_untergrund1,std_b_untergrund1)
print('curvefit',np.sqrt(np.diag(cov_untergrund1)))


x=np.linspace(200,300)
j=0
for i in T1:
    if (gerade(i,*params_untergund1)<0):
        j=j+1
print(j)



plt.figure(1)
plt.plot(T1,I1,'bx',label=r'Messwerte')
plt.plot(T1[:25],(I1[:25]),'rx',label=r'Messpunkte für lineare Reg.')
plt.plot(T1[55:],(I1[55:]),'gx',label=r'Messpunkte für Untergrund')
#plt.plot(x,gerade(x,*params_untergund1),'c-',label=r'Untergrund fit')
#plt.plot(x,gerade(x,m_untergrund1,b_untergrund1),'k-',label=r'Untergrund fit linregress')
# plt.plot(T1[:j],I1[:j],'mx',label=r'Messpunkte ohne Untergrund')
# plt.plot(T1[j:],I1[j:]-gerade(T1[j:],*params_untergund1),'mx',)
# plt.plot(T1,I1,'kx',label=r'Messpunkte mit Untergrund')
plt.legend(loc='best')
plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot1_messwerte.pdf')

I1_richtig=np.append(I1[:j],I1[j:]-gerade(T1[j:],*params_untergund1))


plt.figure(2)
#plt.plot(T1,I1,'bx',label=r'Graph')
plt.plot(x,gerade(x,*params_untergund1),'c-',label=r'Untergrund fit')
#plt.plot(x,gerade(x,m_untergrund1,b_untergrund1),'k-',label=r'Untergrund fit linregress')
plt.plot(T1,I1,'kx',label=r'Messpunkte mit Untergrund')
plt.plot(T1,I1_richtig,'mx',label=r'Messpunkte ohne Untergrund')
plt.legend(loc='best')
plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot1_korrektur.pdf')


#Plot 2


params_untergund2 , cov_untergrund2 = curve_fit(gerade ,T2[55:] , I2[55:])

j=0
for i in T2:
    if (gerade(i,*params_untergund2)<0):
        j=j+1
print(j)



I2_richtig=np.append(I2[:j],I2[j:]-gerade(T2[j:],*params_untergund2))


plt.figure(3)
plt.plot(T2,I2,'bx',label=r'Messwerte')
plt.plot(T2[:25],(I2[:25]),'rx',label=r'Gerade')
plt.plot(T2[53:],(I2[53:]),'gx',label=r'Untergrund')
# plt.plot(x,gerade(x,*params_untergund2),'c-',label=r'Untergrund fit')
# plt.plot(T2[:j],I2[:j],'mx',label=r'Messpunkte ohne Untergrund')
# plt.plot(T2[j:],I2[j:]-gerade(T2[j:],*params_untergund2),'mx',)
plt.legend(loc='best')
plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot2_messwerte.pdf')

plt.figure(4)
plt.plot(T2,I2,'bx',label=r'Messwerte')
# plt.plot(T2[:25],(I2[:25]),'rx',label=r'Messpunkte für lineare Reg.')
# plt.plot(T2[53:],(I2[53:]),'gx',label=r'Messpunkte für Untergrund')
plt.plot(x,gerade(x,*params_untergund2),'c-',label=r'Untergrund fit')
plt.plot(T2,I2_richtig,'mx',label=r'Messpunkte ohne Untergrund')
plt.legend(loc='best')
plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot2_korrektur.pdf')



params_fit1 , cov_fit1 = curve_fit(gerade ,1/T1[:20] ,np.log((I1_richtig[:20])))


plt.figure(5)
plt.plot(1/(T1[:20]),np.log((I1_richtig[:20])),'mx',label=r'log Graph')
plt.plot(1/(T1[:20]),gerade(1/(T1[:20]),*params_fit1),'c-')
#plt.plot(1/(T2[:20]),np.log((I2[:20])),'bx',label=r'log Graph')

plt.legend(loc='best')
plt.xlabel(r' $1/T  1 / \si{\kelvin}$')
plt.ylabel(r'$\log(I)$')
plt.savefig('build/plot3_messung1.pdf')



params_fit2 , cov_fit2 = curve_fit(gerade ,1/T2[:20] ,np.log((I2_richtig[:20])))


plt.figure(6)
plt.plot(1/(T2[:20]),np.log((I2_richtig[:20])),'mx',label=r'log Graph')
plt.plot(1/(T2[:20]),gerade(1/(T2[:20]),*params_fit2),'c-')
#plt.plot(1/(T2[:20]),np.log((I2[:20])),'bx',label=r'log Graph')
plt.legend(loc='best')
plt.xlabel(r' $1/T  1 / \si{\kelvin}$')
plt.ylabel(r'$\log(I)$')
plt.savefig('build/plot3_messung2.pdf')


print('W bei erster messung',params_fit1[0]*const.k/const.e)
print(np.sqrt(np.diag(cov_fit2)))


print('W bei erster messung',params_fit2[0]*const.k/const.e)
print(np.sqrt(np.diag(cov_fit2)))



plt.figure(7)
#plt.plot(T1,I1,'bx',label=r'Graph')
#plt.plot(x,gerade(x,*params_untergund1),'c-',label=r'Untergrund fit')
#plt.plot(x,gerade(x,m_untergrund1,b_untergrund1),'k-',label=r'Untergrund fit linregress')
#plt.plot(T1,I1,'kx',label=r'Messpunkte mit Untergrund')
plt.fill_between(T1,0,I1_richtig,color='m',label=r'Fläche unter den Messpunkten',alpha=0.3)
plt.legend(loc='best')
plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot4_int1.pdf')


plt.figure(8)
#plt.plot(T2,I2,'bx',label=r'Messwerte')
# plt.plot(T2[:25],(I2[:25]),'rx',label=r'Messpunkte für lineare Reg.')
# plt.plot(T2[53:],(I2[53:]),'gx',label=r'Messpunkte für Untergrund')
#plt.plot(x,gerade(x,*params_untergund2),'c-',label=r'Untergrund fit')
#plt.plot(T2,I2_richtig,'mx',label=r'Messpunkte ohne Untergrund')

plt.fill_between(T2,0,I2_richtig,color='c', label=r'Fläche unter den Messpunkten',alpha=0.3)

plt.tight_layout()

plt.legend(loc='best')
plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot4_int2.pdf')
