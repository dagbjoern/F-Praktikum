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



def test(x,a,b,c,d):
    return a*np.abs(np.sin(b*x+c))+d


def kontrast(I_max,I_min):
        k=(I_max-I_min)/(I_min+I_max)
        return k

def brechgas(M,lam,L):
    n=M*lam/L+1
    return n

def brechglas1(M,theta,lam,T):
    test=10/360 * 2*np.pi
    cool=(test+theta)**2-(test-theta)**2
    return 1/(1-M*lam/(T*cool))

def brechglas2(M,theta,lam,T):
    test=10/360 * 2*np.pi
    theata=np.sqrt((test+theta)**2-(test-theta)**2)
    a=M*lam/(2*T)
    return (a**2+2*(1-np.cos(theta))*(1-a))/(2*(1-np.cos(theta)-a))

plt.figure(1)
plt.plot(1,1)
plt.savefig('build/plot1.pdf')

#kontrastwert:
grad , I_max , I_min = np.genfromtxt('Messwertewinkel.txt',unpack=True)
M_2 , M_4 , M_6 , M_8 =np.genfromtxt('messungglas.txt',unpack=True)


M_2mittel=unp.uarray(np.mean(M_2),np.std(M_2))
M_4mittel=unp.uarray(np.mean(M_4),np.std(M_4))
M_6mittel=unp.uarray(np.mean(M_6),np.std(M_6))
M_8mittel=unp.uarray(np.mean(M_8),np.std(M_8))


print('mittelwerte')
print(M_2mittel)
print(M_4mittel)
print(M_6mittel)
print(M_8mittel)


M=np.array([M_2,M_4,M_6,M_8])
M_mittel=unp.uarray([noms(M_2mittel),noms(M_4mittel),noms(M_6mittel),noms(M_8mittel)],[stds(M_2mittel),stds(M_4mittel),stds(M_6mittel),stds(M_8mittel)])


print(M_mittel)
print(M)

phi=np.array([2,4,6,8])
phi=phi
print(phi)
phi=(phi/360) * 2*np.pi
print(phi)
off=10/360 *2*np.pi


messgas=np.array([41,42,42])
# I_min=np.min(I)
# I_max=np.max(I)

winkel=grad/360*2*np.pi

#print('kontrast k',kontrast(I_max,I_min))



#print(tabulate({"Grad": grad,"Winkel": np.round(winkel, 3)," I_min":I_min,"I_max":I_max,"kontrast": np.round(kontrast(I_max,I_min),3) },headers="keys",tablefmt="latex"))


params , cov = curve_fit(test ,winkel,kontrast(I_max,I_min),p0=[0.9,2*np.pi/3,0,0])



print('params',params)
print('fehler',np.sqrt(np.diag(cov)))

param=unp.uarray(params,np.sqrt(np.diag(cov)))
print((np.pi/2-param[2])/param[1])

print(((np.pi/2-param[2])/param[1])/(2*np.pi)*360)




x=np.linspace(-1,4,1000)

plt.figure(1)
plt.plot(winkel,kontrast(I_max,I_min),'kx',label=r'Messwerte')
plt.plot(x,test(x,*params),label=r'Fit-fkt.')
plt.legend(loc='best')
plt.xlabel(r'Winkel $\Phi/ rad$')
plt.ylabel(r'Kontrast $K$')
plt.savefig('build/plot1.pdf')

##############################################
print(brechgas(messgas,632.990*10**(-9),0.1))
print(np.mean(brechgas(messgas,632.990*10**(-9),0.1)),np.std(brechgas(messgas,632.990*10**(-9),0.1)))

print('Abweichung',(np.mean(brechgas(messgas,632.990*10**(-9),0.1))-1.000292)/1.000292)
print(np.size(M[1,:]))

plt.figure(2)
plt.errorbar(phi,noms(M_mittel), xerr=0, yerr=std(M_mittel) ,'kx',label=r'Messwerte')
plt.legend(loc='best')
plt.xlabel(r'Winkel $\Phi/ rad$')
plt.ylabel(r'M $K$')
plt.savefig('build/plotglas.pdf')




#def brechglas(M,theta,lam,T):

print(brechglas1(M_2,phi[0],632.990*10**(-9),0.001))
print(brechglas2(M_2,phi[0],632.990*10**(-9),0.001))
print(brechglas1(M_4,phi[1],632.990*10**(-9),0.001))
print(brechglas2(M_4,phi[1],632.990*10**(-9),0.001))
print(brechglas1(M_6,phi[2],632.990*10**(-9),0.001))
print(brechglas2(M_6,phi[2],632.990*10**(-9),0.001))
print(brechglas1(M_8,phi[3],632.990*10**(-9),0.001))
print(brechglas2(M_8,phi[3],632.990*10**(-9),0.001))
#
# print(tabulate({"M_2": M_2 ,"n_2": np.round(brechglas2(M_2,phi[0],632.990*10**(-9),0.001), 2),"M_4":M_4,"n_4":np.round(brechglas2(M_4,phi[1],632.990*10**(-9),0.001), 2),"M_6":M_6,"n_6":np.round(brechglas2(M_6,phi[2],632.990*10**(-9),0.001), 2),"M_8":M_8,"n_8":np.round(brechglas2(M_8,phi[3],632.990*10**(-9),0.001), 2) },headers="keys",tablefmt="latex"))
#
#
# print(np.mean([brechglas2(M_2,phi[0],632.990*10**(-9),0.001),brechglas2(M_4,phi[1],632.990*10**(-9),0.001),brechglas2(M_6,phi[2],632.990*10**(-9),0.001),brechglas2(M_8,phi[3],632.990*10**(-9),0.001)]))
# print(np.std([brechglas2(M_2,phi[0],632.990*10**(-9),0.001),brechglas2(M_4,phi[1],632.990*10**(-9),0.001),brechglas2(M_6,phi[2],632.990*10**(-9),0.001),brechglas2(M_8,phi[3],632.990*10**(-9),0.001)]))

print('Abweichung glas',(np.mean([brechglas2(M_2,phi[0],632.990*10**(-9),0.001),brechglas2(M_4,phi[1],632.990*10**(-9),0.001),brechglas2(M_6,phi[2],632.990*10**(-9),0.001),brechglas2(M_8,phi[3],632.990*10**(-9),0.001)])-1.45 )/1.45 )
