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


#kontrastwert:
grad , I_max , I_min = np.genfromtxt('Messwertewinkel.txt',unpack=True)
#
# I_min=np.min(I)
# I_max=np.max(I)

winkel=grad/360*2*np.pi

#print('kontrast k',kontrast(I_max,I_min))



print(tabulate({"Grad": grad,"Winkel": np.round(winkel, 3)," I_min":I_min,"I_max":I_max,"kontrast": np.round(kontrast(I_max,I_min),3) },headers="keys",tablefmt="latex"))


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
# plt.xlabel(r'Temperatur $T/ \si{\kelvin}$')
# plt.ylabel(r'Strom $I/ \si{\ampere} $')
plt.savefig('build/plot1.pdf')

##############################################
