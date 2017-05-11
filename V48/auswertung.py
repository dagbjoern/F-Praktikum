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
print('params',params_untergund1)
print('fehler',np.sqrt(np.diag(cov_untergrund1)))





x=np.linspace(200,300)
j=0
for i in T1:
    if (gerade(i,*params_untergund1)<0):
        j=j+1
print(j)



plt.figure(1)
plt.plot(T1,I1,'kx',label=r'Messwerte')
#plt.plot(T1[:25],(I1[:25]),'rx',label=r'Fit')
plt.plot(T1[55:],(I1[55:]),'bx',label=r'Untergrund')
plt.plot(T1[:40],(I1[:40]),'bx',label=r'test')

#plt.plot(x,gerade(x,*params_untergund1),'c-',label=r'Untergrund fit')
#plt.plot(x,gerade(x,m_untergrund1,b_untergrund1),'k-',label=r'Untergrund fit linregress')
# plt.plot(T1[:j],I1[:j],'mx',label=r'Messpunkte ohne Untergrund')
# plt.plot(T1[j:],I1[j:]-gerade(T1[j:],*params_untergund1),'mx',)
# plt.plot(T1,I1,'kx',label=r'Messpunkte mit Untergrund')
plt.legend(loc='best')
#plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
#plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot1_messwerte.pdf')

I1_richtig=np.append(I1[:j],I1[j:]-gerade(T1[j:],*params_untergund1))


plt.figure(2)
#plt.plot(T1,I1,'bx',label=r'Graph')
plt.plot(x,gerade(x,*params_untergund1),'c-',label=r'Untergrund-Fit')
#plt.plot(x,gerade(x,m_untergrund1,b_untergrund1),'k-',label=r'Untergrund fit linregress')
plt.plot(T1,I1,'kx',label=r'Messpunkte')
plt.plot(T1,I1_richtig,'mx',label=r'ohne Untergrund')
plt.legend(loc='best')
#plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
#plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot1_korrektur.pdf')


#Plot 2


params_untergund2 , cov_untergrund2 = curve_fit(gerade ,T2[55:] , I2[55:])

print('params',params_untergund2)
print('fehler',np.sqrt(np.diag(cov_untergrund2)))



j=0
for i in T2:
    if (gerade(i,*params_untergund2)<0):
        j=j+1
print(j)



I2_richtig=np.append(I2[:j],I2[j:]-gerade(T2[j:],*params_untergund2))


plt.figure(3)
plt.plot(T2,I2,'kx',label=r'Messwerte')
#plt.plot(T2[:25],(I2[:25]),'rx',label=r'Fit')
plt.plot(T2[53:],(I2[53:]),'bx',label=r'Untergrund')
# plt.plot(x,gerade(x,*params_untergund2),'c-',label=r'Untergrund fit')
# plt.plot(T2[:j],I2[:j],'mx',label=r'Messpunkte ohne Untergrund')
# plt.plot(T2[j:],I2[j:]-gerade(T2[j:],*params_untergund2),'mx',)
plt.legend(loc='best')
#plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
#plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot2_messwerte.pdf')

plt.figure(4)
# plt.plot(T2[:25],(I2[:25]),'rx',label=r'Messpunkte für lineare Reg.')
# plt.plot(T2[53:],(I2[53:]),'gx',label=r'Messpunkte für Untergrund')
plt.plot(x,gerade(x,*params_untergund2),'c-',label=r'Untergrund-Fit')
plt.plot(T2,I2,'kx',label=r'Messwerte')
plt.plot(T2,I2_richtig,'mx',label=r'ohne Untergrund')
plt.legend(loc='best')
#plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
#plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot2_korrektur.pdf')



params_fit1 , cov_fit1 = curve_fit(gerade ,1/T1[:20] ,np.log((I1_richtig[:20])))


print('params cool1',params_fit1)
print('fehler cool1',np.sqrt(np.diag(cov_fit1)))
m_1=unp.uarray(params_fit1[0],np.sqrt(np.diag(cov_fit1)[0]))
b_1=unp.uarray(params_fit1[1],np.sqrt(np.diag(cov_fit1)[0]))
print('m_1',m_1,' \n ','b_1',b_1)



plt.figure(5)
plt.plot(1/(T1[:20]),np.log((I1_richtig[:20])),'gx',label=r'Näherung')
plt.plot(1/(T1[:20]),gerade(1/(T1[:20]),*params_fit1),'c-',label=r'Fit')
#plt.plot(1/(T2[:20]),np.log((I2[:20])),'bx',label=r'log Graph')

plt.legend(loc='best')
#plt.xlabel(r' $1/T  1 / \si{\kelvin}$')
#plt.ylabel(r'$\log(I)$')
plt.savefig('build/plot3_messung1.pdf')



params_fit2 , cov_fit2 = curve_fit(gerade ,1/T2[:20] ,np.log((I2_richtig[:20])))


print('params cool2',params_fit2)
print('fehler cool2',np.sqrt(np.diag(cov_fit2)))
m_2=unp.uarray(params_fit2[0],np.sqrt(np.diag(cov_fit2)[0]))
b_2=unp.uarray(params_fit2[1],np.sqrt(np.diag(cov_fit2)[0]))
print('m_1',m_2,'/n','b_2',b_2)


plt.figure(6)
plt.plot(1/(T2[:20]),np.log((I2_richtig[:20])),'gx',label=r'Näherung')
plt.plot(1/(T2[:20]),gerade(1/(T2[:20]),*params_fit2),'c-',label=r'Fit')
#plt.plot(1/(T2[:20]),np.log((I2[:20])),'bx',label=r'log Graph')
plt.legend(loc='best')
#plt.xlabel(r' $1/T  1 / \si{\kelvin}$')
#plt.ylabel(r'$\log(I)$')
plt.savefig('build/plot3_messung2.pdf')




print('W bei erster messung',-m_1*const.k/const.e)


print('W bei zweiten messung',-m_2*const.k/const.e)


plt.figure(7)
plt.plot(T1,I1_richtig,'mx',label=r'Korr. Messwerte')
plt.plot(T1[:25],(I1_richtig[:25]),'gx',label=r'Näherung')
#plt.plot(T2[53:],(I2[53:]),'bx',label=r'Untergrund')
# plt.plot(x,gerade(x,*params_untergund2),'c-',label=r'Untergrund fit')
# plt.plot(T2[:j],I2[:j],'mx',label=r'Messpunkte ohne Untergrund')
# plt.plot(T2[j:],I2[j:]-gerade(T2[j:],*params_untergund2),'mx',)
plt.legend(loc='best')
#plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
#plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot1_fit.pdf')




plt.figure(8)
plt.plot(T2,I2_richtig,'mx',label=r'Korr. Messwerte')
plt.plot(T2[:25],(I2_richtig[:25]),'gx',label=r'Näherung')
#plt.plot(T2[53:],(I2[53:]),'bx',label=r'Untergrund')
# plt.plot(x,gerade(x,*params_untergund2),'c-',label=r'Untergrund fit')
# plt.plot(T2[:j],I2[:j],'mx',label=r'Messpunkte ohne Untergrund')
# plt.plot(T2[j:],I2[j:]-gerade(T2[j:],*params_untergund2),'mx',)
plt.legend(loc='best')
#plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
#plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot2_fit.pdf')





plt.figure(9)
#plt.plot(T1,I1,'bx',label=r'Graph')
#plt.plot(x,gerade(x,*params_untergund1),'c-',label=r'Untergrund fit')
#plt.plot(x,gerade(x,m_untergrund1,b_untergrund1),'k-',label=r'Untergrund fit linregress')
#plt.plot(T1,I1,'kx',label=r'Messpunkte mit Untergrund')
plt.fill_between(T1,0,I1_richtig,label=r'Fläche',alpha=0.4,hatch='//',edgecolor='k')
plt.legend(loc='best')
#plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
#plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot4_int1.pdf')


plt.figure(10)
#plt.plot(T2,I2,'bx',label=r'Messwerte')
# plt.plot(T2[:25],(I2[:25]),'rx',label=r'Messpunkte für lineare Reg.')
# plt.plot(T2[53:],(I2[53:]),'gx',label=r'Messpunkte für Untergrund')
#plt.plot(x,gerade(x,*params_untergund2),'c-',label=r'Untergrund fit')
#plt.plot(T2,I2_richtig,'mx',label=r'Messpunkte ohne Untergrund')

plt.fill_between(T2,0,I2_richtig,label=r'Fläche',alpha=0.4,hatch='//',edgecolor='k')

plt.tight_layout()

plt.legend(loc='best')
#plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
#plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot4_int2.pdf')


A1= trapz(I1_richtig,T1)
A2= trapz(I2_richtig,T2)

#heizraten
b_h1=1.2/60
b_h2=2/60

print('Flache 1', A1)
print('Flache 2', A2)


A_1=np.array(I1_richtig)

for index, value in enumerate(I1_richtig):
    A_1[index]= trapz(I1_richtig[index:],T1[index:])

print(A_1)




params_reg1 , cov_reg1 = curve_fit(gerade ,1/T1[:55],np.log(A_1[:55]/(I1_richtig[:55])))

m_1reg=unp.uarray(params_reg1[0],np.sqrt(np.diag(cov_reg1)[0]))
b_1reg=unp.uarray(params_reg1[1],np.sqrt(np.diag(cov_reg1)[0]))
print('m_1reg',m_1reg,'/n','b_1reg',b_1reg)
print('W bei erster messung mit intergral',m_1reg*const.k/const.e)


plt.figure(11)
plt.plot(1/T1[:len(A_1)-1],np.log(A_1[:len(A_1)-1]/(I1_richtig[:len(A_1)-1])),'kx',label=r'Messwerte')
plt.plot(1/T1[:55],np.log(A_1[:55]/(I1_richtig[:55])),'rx',label=r'benutze Messwerte')

# plt.plot(T2[53:],(I2[53:]),'gx',label=r'Messpunkte für Untergrund')
plt.plot(1/T1[:55],gerade(1/T1[:55],*params_reg1),'c-',label=r'Fit')
#plt.plot(T2,I2_richtig,'mx',label=r'Messpunkte ohne Untergrund')
plt.legend(loc='best')
#plt.xlabel(r'$\frac{1}{T} \si{\per\kelvin}$')
#plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot5_1.pdf')




A_2=np.array(I2_richtig)

for index, value in enumerate(I2_richtig):
    A_2[index]= trapz(I2_richtig[index:],T2[index:])



params_reg2 , cov_reg2 = curve_fit(gerade ,1/T2[:55],np.log(np.abs(A_2[:55]/I2_richtig[:55])))

m_2reg=unp.uarray(params_reg2[0],np.sqrt(np.diag(cov_reg2)[0]))
b_2reg=unp.uarray(params_reg2[1],np.sqrt(np.diag(cov_reg2)[0]))
print('m_2',m_2reg,'/n','b_2',b_2reg)

print('W bei zweiter messung mit intergral',m_2reg*const.k/const.e)



plt.figure(12)
plt.plot(1/T2[:len(A_2)-1],np.log(A_2[:len(A_2)-1]/(I2_richtig[:len(A_2)-1])),'kx',label=r'Messwerte')
plt.plot(1/T2[:55],np.log(A_2[:55]/(I2_richtig[:55])),'rx',label=r'Benutze Messwerte')
# plt.plot(T2[53:],(I2[53:]),'gx',label=r'Messpunkte für Untergrund')
plt.plot(1/T2[:len(A_2)-1],gerade(1/T2[:len(A_2)-1],*params_reg2),'c-',label=r'Fit')
#plt.plot(T2,I2_richtig,'mx',label=r'Messpunkte ohne Untergrund')
plt.legend(loc='best')
#plt.xlabel(r'$\frac{1}{T} \si{\per\kelvin}$')
#plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot5_2.pdf')








#print(tabulate({"Strom":I1_richtig,"Temperatur":T1},tablefmt="latex"))
