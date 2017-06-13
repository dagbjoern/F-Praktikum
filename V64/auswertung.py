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
    return 1/(1-M*lam/(T*theta**2))

def brechglas2(M,theta,lam,T):
        a=M*lam/(2*T)
        return (a**2+2*(1-np.cos(theta))*(1-a))/(2*(1-np.cos(theta)-a))

#kontrastwert:
grad , I_max , I_min = np.genfromtxt('Messwertewinkel.txt',unpack=True)
M_2 , M_4 , M_6 , M_8 =np.genfromtxt('messungglas.txt',unpack=True)


phi=np.array([2,4,6,8])
phi=phi+10
print(phi)
phi=(phi/360) * 2*np.pi
print(phi)

messgas=np.array([41,42,42])
# I_min=np.min(I)
# I_max=np.max(I)

winkel=grad/360*2*np.pi

#print('kontrast k',kontrast(I_max,I_min))



#print(tabulate({"Grad": grad,"Winkel": np.round(winkel, 3)," I_min":I_min,"I_max":I_max,"kontrast": np.round(kontrast(I_max,I_min),3) },headers="keys",tablefmt="latex"))


params , cov = curve_fit(test ,winkel,kontrast(I_max,I_min),p0=[0.9,2*np.pi/3,0,0])



print('params',params)
print('fehler',np.sqrt(np.diag(cov)))

print((np.pi/2-params[2])/params[1])

print(((np.pi/2-params[2])/params[1])/(2*np.pi)*360)




x=np.linspace(-1,4,1000)

plt.figure(1)
plt.plot(winkel,kontrast(I_max,I_min),'kx',label=r'Messwerte')
plt.plot(x,test(x,*params))
plt.legend(loc='best')
#plt.xlabel(r'Winkel $\Phi/ \si{\radian}$')
#plt.ylabel(r'Kontrast $K$')
plt.savefig('build/plot1.pdf')

##############################################
print(brechgas(messgas,632.990*10**(-9),0.1))
print(np.mean(brechgas(messgas,632.990*10**(-9),0.1)),np.std(brechgas(messgas,632.990*10**(-9),0.1)))

#def brechglas(M,theta,lam,T):
print('glas')
#print(brechglas1(M_2,phi[0],632.990*10**(-9),0.001))
print(brechglas2(M_2,phi[0],632.990*10**(-9),0.001))
#print(brechglas1(M_4,phi[1],632.990*10**(-9),0.001))
print(brechglas2(M_4,phi[1],632.990*10**(-9),0.001))
#print(brechglas1(M_6,phi[2],632.990*10**(-9),0.001))
print(brechglas2(M_6,phi[2],632.990*10**(-9),0.001))
#print(brechglas1(M_8,phi[3],632.990*10**(-9),0.001))
print(brechglas2(M_8,phi[3],632.990*10**(-9),0.001))

#print(tabulate({"Grad": grad,"Winkel": np.round(winkel, 3)," I_min":I_min,"I_max":I_max,"kontrast": np.round(kontrast(I_max,I_min),3) },headers="keys",tablefmt="latex"))
