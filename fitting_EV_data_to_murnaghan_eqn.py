from pylab import *
from scipy.optimize import leastsq
import io
import os

# Figure attributes
plt.rcParams.update({'font.size': 14})
plt.rcParams['figure.figsize'] = 12, 10

# Reading data from file 'Output.txt'
path='./Output.txt'
if os.path.exists(path):
  bcc_data=np.loadtxt('Output.txt')
else:
  print('The data file is missing!')

#bcc_data=np.loadtxt('Output.txt')
volumes=np.divide(bcc_data[:,1],2)
energies=np.divide(bcc_data[:,0],2)
fit_volumes = np.linspace(min(volumes),max(volumes),100)

p,q,r=polyfit(volumes,energies,2) 
#Initial guesses
v_ini = -q/(2*p)
e_ini = (p*v_ini**2) + (q*v_ini + r)
b_ini = 2*p*v_ini
b_prime = 4

# List made up of initial guesses 
initial_guesses = [e_ini, b_ini, b_prime, v_ini]

''' Murnaghan equation is taken from: First-principles calculation of the equilibrium ground-state properties of transition metals: Applications to Nb and Mo authored by CL Fu, KM Ho - Physical Review B, 1983 '''  
def murnaghan_eqn(parameters,volume):
    E_ini, B_ini, B_prime, V_ini = parameters[0], parameters[1], parameters[2], parameters[3]
    Energy = E_ini + B_ini*volume/B_prime*(((V_ini/volume)**B_prime)/(B_prime-1.0)+1.0) - V_ini*B_ini/(B_prime-1.0)
    return Energy

# Error function:
def error_func(pars,y,x):
    error =  y - murnaghan_eqn(pars,x)       #error function to be minimized
    return error

murnaghan_fit_pars, ier = leastsq(error_func, initial_guesses, args=(energies,volumes)) 

#Plotting:
plot(volumes,energies,'bx', label = 'Interatomic potential data points')
plot(fit_volumes, murnaghan_eqn(murnaghan_fit_pars,fit_volumes), 'g', label='Murnaghan fit')

ax = gca()
text(0.2,0.5,'Minimum volume = %1.2f $\AA^3$' % murnaghan_fit_pars[3], transform = ax.transAxes, color = 'magenta')
text(0.2,0.3,'Minimum energy = %1.2f eV/atom' % murnaghan_fit_pars[0], transform = ax.transAxes, color = 'magenta')
text(0.2,0.4,'Bulk modulus = %1.2f GPa' % (murnaghan_fit_pars[1]*160.21773), transform = ax.transAxes, color = 'magenta')

xlabel('Atomic volume ($\AA^3$)')
ylabel('Atomic energy (eV/atom)')
title('Fe energy-volume data fit to Murnaghan equation',loc='center',color='red',fontsize=24)
legend(loc='best',fontsize=18)
savefig('Fe_bcc_bulk.png')
show()


