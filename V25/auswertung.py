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


xB_0, IB_0 = np.genfromtxt('OhneBFeld.txt',unpack=True)
x1, I1 = np.genfromtxt('0.4Ampere.txt',unpack=True)
x2, I2 = np.genfromtxt('0.5Ampere.txt',unpack=True)
x3, I3 = np.genfromtxt('0.6Ampere.txt',unpack=True)
x4, I4 = np.genfromtxt('0.735Ampere.txt',unpack=True)
x5, I5 = np.genfromtxt('0.8Ampere.txt',unpack=True)
x6, I6 = np.genfromtxt('0.91Ampere.txt',unpack=True)
#x7, I7 = np.genfromtxt('0.4Ampere.txt',unpack=True)





#
# def fehler_b(std_m,x):
#     return(std_m*np.sqrt(np.mean(x**2)))


def gerade(x,m,b):
    return m*x+b

#linerare regression

# m_untergrund1 , b_untergrund1 , r ,p ,std_m_untergrund1 =stats.linregress(T1[55:],(I1[55:]))

# std_b_untergrund1=fehler_b(std_m_untergrund1,T1[55:])



# print('lin ',std_m_untergrund1,std_b_untergrund1)
#print('params',params_untergund1)
#print('fehler',np.sqrt(np.diag(cov_untergrund1)))





#x=np.linspace(200,300)
#j=0
#for i in T1:
#    if (gerade(i,*params_untergund1)<0):
#        j=j+1
#print(j)
def gauss(x,a,b,c,d):
    return a*np.exp(-(x-b)**2/(2*c**2))+d

params_B0 , cov_0 = curve_fit(gauss,xB_0[40:52],IB_0[40:52],p0=[850,7,0.5,150])

fehler_0=np.sqrt(np.diag(cov_0))
params_0 = unp.uarray(params_B0,np.sqrt(np.diag(cov_0)))

print('Nullpositon=',params_0[1])

plt.figure(1)
plt.plot(xB_0, IB_0,'kx',alpha=0.5,label=r'Messwerte')
plt.plot(xB_0[40:52],gauss(xB_0[40:52],*params_B0),'-',label=r'Fit')
#plt.plot(xB_0,gauss(xB_0,900,7,0.5),'--r',label=r'Fit')

plt.legend(loc='best')
plt.xlabel(r'Spindelumdrehungen $Skt$')
plt.ylabel(r'Strom $I/ \si{\micro\volt} $')
plt.tight_layout()
plt.savefig('build/plot_B=0.pdf')

a=0.0025  #radius der Polschuhe
L=0.070   #Länge der Polschuhe
Skt=0.0018#eine spindelumdrehung
Abstand=np.zeros([6])
Strom=np.array([0.4,0.5,0.6,0.735,0.8,0.91])
B_Feld=np.array([0.3,0.38,0.45,0.55,0.595,0.765])
delta_B=0.968*B_Feld/a
print('Delta_B',delta_B)
#params_1 , cov_1 = curve_fit(gerade ,(T1[55:]),(I1[55:]))
params_1u , cov_1u = curve_fit(gauss,x1[30:46],I1[30:46],p0=[400,6,0.5,150])
params_1o , cov_1o = curve_fit(gauss,x1[54:68],I1[54:68],p0=[400,8,0.5,150])
params_1unten = unp.uarray(params_1u,np.sqrt(np.diag(cov_1u)))
params_1oben = unp.uarray(params_1o,np.sqrt(np.diag(cov_1o)))
Abstand_1 = params_1oben[1] - params_1unten[1]
s_1=np.array([params_1oben[1]-params_0[1],params_0[1] - params_1unten[1],Abstand_1/2])
s_1=s_1
print('b_l 1',params_1unten[1])
print('b_r 1',params_1oben[1])
print('abstand O  1',params_1oben[1]-params_0[1])
print('abstand U  1',params_0[1] - params_1unten[1])
print('abstand /2 1',Abstand_1/2)

s_1_mean=unp.uarray([np.mean(noms(s_1))],[np.std(noms(s_1))])



plt.figure(2)
plt.plot(x1,I1,'kx',alpha=0.5,label=r'Messwerte')
plt.plot(x1[30:47],gauss(x1[30:47],*params_1u),'-',label=r'Fit 1')
plt.plot(x1[54:68],gauss(x1[54:68],*params_1o),'-',label=r'Fit 2')
#plt.plot(x,gerade(x,m_untergrund1,b_untergrund1),'k-',label=r'Untergrund fit linregress')
plt.legend(loc='best')
plt.xlabel(r'Spindelumdrehungen $Skt$')
plt.ylabel(r'Strom $I/ \si{\micro\volt} $')
plt.tight_layout()
plt.savefig('build/plot1.pdf')


params_2u , cov_2u = curve_fit(gauss,x2[27:43],I2[27:43],p0=[400,6,0.5,150])
params_2o , cov_2o = curve_fit(gauss,x2[55:70],I2[55:70],p0=[400,8,0.5,150])

params_2unten = unp.uarray(params_2u,np.sqrt(np.diag(cov_2u)))
params_2oben = unp.uarray(params_2o,np.sqrt(np.diag(cov_2o)))
Abstand_2 = params_2oben[1] - params_2unten[1]
s_2=np.array([params_2oben[1]-params_0[1],params_0[1] - params_2unten[1],Abstand_2/2])
s_2=s_2
s_2_mean=unp.uarray(np.mean(noms(s_2)),np.std(noms(s_2)))
print('b_l 2',params_2unten[1])
print('b_r 2',params_2oben[1])

print('abstand O 2',params_2oben[1]-params_0[1])
print('abstand U 2',params_0[1] - params_2unten[1])
print('abstand /2  2',Abstand_2/2)




plt.figure(3)
plt.plot(x2, I2,'kx',alpha=0.5,label=r'Messwerte')
plt.plot(x2[27:43],gauss(x2[27:43],*params_2u),'-',label=r'Fit 1')
plt.plot(x2[55:70],gauss(x2[55:70],*params_2o),'-',label=r'Fit 2')
plt.legend(loc='best')
plt.xlabel(r'Spindelumdrehungen $Skt$')
plt.ylabel(r'Strom $I/ \si{\micro\volt} $')
plt.tight_layout()
plt.savefig('build/plot2.pdf')

params_3u , cov_3u = curve_fit(gauss,x3[33:45],I3[33:45],p0=[350,6,0.5,10])
params_3o , cov_3o = curve_fit(gauss,x3[62:75],I3[62:75],p0=[350,7.5,0.5,10])

params_3unten = unp.uarray(params_3u,np.sqrt(np.diag(cov_3u)))
params_3oben = unp.uarray(params_3o,np.sqrt(np.diag(cov_3o)))
Abstand_3 = params_3oben[1] - params_3unten[1]
s_3=np.array([params_3oben[1]-params_0[1],params_0[1] - params_3unten[1],Abstand_3/2])
s_3=s_3
s_3_mean=unp.uarray([np.mean(noms(s_3))],[np.std(noms(s_3))])
print('b_l 3',params_3unten[1])
print('b_r 3',params_3oben[1])

print('abstand O 3',params_3oben[1]-params_0[1])
print('abstand U 3',params_0[1] - params_3unten[1])
print('abstand /2  3',Abstand_3/2)

plt.figure(4)
plt.plot(x3, I3,'kx',alpha=0.5,label=r'Messwerte')
plt.plot(x3[33:45],gauss(x3[33:45],*params_3u),'-',label=r'Fit 1')
plt.plot(x3[62:75],gauss(x3[62:75],*params_3o),'-',label=r'Fit 2')
plt.legend(loc='best')
plt.xlabel(r'Spindelumdrehungen $Skt$')
plt.ylabel(r'Strom $I/ \si{\micro\volt} $')
plt.tight_layout()
plt.savefig('build/plot3.pdf')




params_4u , cov_4u = curve_fit(gauss,x4[29:45],I4[29:45],p0=[350,6,0.5,10])
params_4o , cov_4o = curve_fit(gauss,x4[62:80],I4[62:80],p0=[350,7.5,0.5,10])

params_4unten = unp.uarray(params_4u,np.sqrt(np.diag(cov_4u)))
params_4oben = unp.uarray(params_4o,np.sqrt(np.diag(cov_4o)))
Abstand_4 = params_4oben[1] - params_4unten[1]
s_4=np.array([params_4oben[1]-params_0[1],params_0[1] - params_4unten[1],Abstand_4/2])
s_4=s_4
s_4_mean=unp.uarray(np.mean(noms(s_4)),np.std(noms(s_4)))
print('b_l 4',params_4unten[1])
print('b_r 4',params_4oben[1])


print('abstand O 4',params_4oben[1]-params_0[1])
print('abstand U 4',params_0[1] - params_4unten[1])
print('abstand /2  4',Abstand_4/2)


plt.figure(5)
plt.plot(x4, I4,'kx',alpha=0.5,label=r'Messwerte')
plt.plot(x4[29:45],gauss(x4[29:45],*params_4u),'-',label=r'Fit 1')
plt.plot(x4[62:80],gauss(x4[62:80],*params_4o),'-',label=r'Fit 2')
plt.legend(loc='best')
plt.xlabel(r'Spindelumdrehungen $Skt$')
plt.ylabel(r'Strom $I/ \si{\micro\volt} $')
plt.tight_layout()
plt.savefig('build/plot4.pdf')



params_5u , cov_5u = curve_fit(gauss,x5[29:45],I5[29:45],p0=[350,6,0.5,10])
params_5o , cov_5o = curve_fit(gauss,x5[62:78],I5[62:78],p0=[350,7.5,0.5,10])

params_5unten = unp.uarray(params_5u,np.sqrt(np.diag(cov_5u)))
params_5oben = unp.uarray(params_5o,np.sqrt(np.diag(cov_5o)))
Abstand_5 = params_5oben[1] - params_5unten[1]
s_5=np.array([params_5oben[1]-params_0[1],params_0[1] - params_5unten[1],Abstand_5/2])
s_5=s_5
s_5_mean=unp.uarray(np.mean(noms(s_5)),np.std(noms(s_5)))
print('b_l 5',params_5unten[1])
print('b_r 5',params_5oben[1])


print('abstand O 5',params_5oben[1]-params_0[1])
print('abstand U 5',params_0[1] - params_5unten[1])
print('abstand /2  5',Abstand_5/2)


plt.figure(6)
plt.plot(x5, I5,'kx',alpha=0.5,label=r'Messwerte')
plt.legend(loc='best')
plt.plot(x5[29:45],gauss(x5[29:45],*params_5u),'-',label=r'Fit 1')
plt.plot(x5[62:78],gauss(x5[62:78],*params_5o),'-',label=r'Fit 2')
plt.xlabel(r'Spindelumdrehungen $Skt$')
plt.ylabel(r'Strom $I/ \si{\micro\volt} $')
plt.tight_layout()
plt.savefig('build/plot5.pdf')


params_6u , cov_6u = curve_fit(gauss,x6[29:43],I6[29:43],p0=[350,6,0.5,10])
params_6o , cov_6o = curve_fit(gauss,x6[66:83],I6[66:83],p0=[350,7.5,0.5,10])

params_6unten = unp.uarray(params_6u,np.sqrt(np.diag(cov_6u)))
params_6oben = unp.uarray(params_6o,np.sqrt(np.diag(cov_6o)))
Abstand_6 = params_6oben[1] - params_6unten[1]

s_6=np.array([params_6oben[1]-params_0[1],params_0[1] - params_6unten[1],Abstand_6/2])
s_6=s_6
s_6_mean=unp.uarray(np.mean(noms(s_6)),np.std(noms(s_6)))
print('b_l 6',params_6unten[1])
print('b_r 6',params_6oben[1])

print('abstand O 6',params_6oben[1]-params_0[1])
print('abstand U 6',params_0[1] - params_6unten[1])
print('abstand /2  6',Abstand_6/2)



plt.figure(7)
plt.plot(x6, I6,'kx',alpha=0.5,label=r'Messwerte')
plt.plot(x6[29:43],gauss(x6[29:43],*params_6u),'-',label=r'Fit 1')
plt.plot(x6[66:83],gauss(x6[66:83],*params_6o),'-',label=r'Fit 2')
plt.legend(loc='best')
plt.xlabel(r'Spindelumdrehungen $Skt$')
plt.ylabel(r'Strom $I/ \si{\micro\volt} $')
plt.tight_layout()
plt.savefig('build/plot6.pdf')


Abstande=unp.uarray([noms(s_1_mean),noms(s_2_mean),noms(s_3_mean),noms(s_4_mean),noms(s_5_mean),noms(s_6_mean)],[stds(s_1_mean),stds(s_2_mean),stds(s_3_mean),stds(s_4_mean),stds(s_5_mean),stds(s_6_mean)])
print('Abstand mittelwerte',Abstande)

Abstande=Abstande*Skt
print('Abstand mittelwerte in meter',Abstande)

def gerade(x,m):
    T=190+273.15
    L=0.070   #Länge der Polschuhe
    l=0.455   #Länge zu blende
    cool=(l*L*(1-(L)/(2*l)))/(6*const.k*T)
    return m*cool*x


params_lin , cov_lin = curve_fit(gerade,delta_B,noms(Abstande))



plt.figure(8)
plt.errorbar(delta_B, noms(Abstande),xerr=0,yerr=stds(Abstande),LineStyle='none',label=r'Messwerte')
plt.plot(delta_B, gerade(delta_B,*params_lin),label=r'Fit-Funktion')
#plt.plot(delta_B, gerade(delta_B,const.value('Bohr magneton')),label=r'theo')

plt.legend(loc='best')
plt.xlabel(r'$\frac{\partial B}{\partial z}/ \si{\tesla\per\meter}$')
plt.ylabel(r'$\bar{s} / \si{\meter}$')
plt.tight_layout()
plt.savefig('build/lin_plot.pdf')
params_bohr = unp.uarray(params_lin,np.sqrt(np.diag(cov_lin)))

print('mu_B gemessen',params_bohr)
print('mu_B theo', const.value('Bohr magneton'))####
print('Abweichung zur Theorie',(params_bohr-const.value('Bohr magneton'))/const.value('Bohr magneton'))
