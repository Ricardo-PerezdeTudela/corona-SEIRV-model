# zombie apocalypse modeling
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
plt.rcParams['figure.figsize'] = 8, 8

# Constants
L = 271.23 # Influx rate
tcE0 = 3.11e-8 # Transmission between S and E
tcI0 = 0.62e-8 # Transmission between S and I
tcV0 = 1.03e-8 # Transmission between S and V
tac = 1.01e-4 # Transmission adjustment coefficient
#tac=0 #constant transmission rate -> unrealisticly pessimistic model

ndr = 3.01e-5  # natural death rate
ddr = 0.01 # disease-induced death rate

ip = 7.   # Incubation period
rr = 1./15. # Recovery rate

rrv = 1. # Removal rate of virus
srE = 2.30 # shielding rate by exposed people
srI = 0. # shielding rate by infected people


# solve the system dy/dt = f(y, t)
def f(y, t):
    S, E, I, R, V = y

    a = 1./ip
    tcE = tcE0/(1.+tac*E)
    tcI = tcI0/(1.+tac*I)
    tcV = tcV0/(1.+tac*V)
    X = tcE*S*E+tcI*S*I+tcV*S*V

    f0 = L-X-ndr*S
    f1 = X-(a+ndr)*E
    f2 = a*E-(ndr+ddr+rr)*I
    f3 = rr*I-ndr*R
    f4 = srE*E+srI*I-rrv*V

    return [f0, f1, f2, f3, f4]

# initial conditions
S = 8998505
E = 1000
I = 475
R = 10
V = 10000

y0 = [S, E, I, R, V]


t  = np.linspace(0, 180., 181) # time grid


# solve the DEs
soln = odeint(f, y0, t)
S = soln[:, 0]
E = soln[:, 1]
I = soln[:, 2]
R = soln[:, 3]
V = soln[:, 4]

Icum = np.cumsum(I)


plt.figure()
plt.plot(t, I, label='Infected')
plt.plot(t, E, label='Exposed')
#plt.plot(t, I*ddr, label='Dead by the virus')
#plt.plot(t, Icum*ddr, 'bo', label='Dead by the virus (cum)')

plt.plot(t, S*ndr, label='Dead by natural means')


plt.xlabel('Days from outbreak')
plt.ylabel('Population')

plt.xticks(np.arange(min(t), max(t)+1, 30))
plt.grid(color='black', linestyle='--', linewidth=1)

plt.legend(loc=0)
plt.show()


#plt.plot(E, I, label='Phase Space')
#plt.show()

