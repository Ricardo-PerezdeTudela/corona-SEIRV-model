import numpy as np
import matplotlib.pyplot as plt
from gekko import GEKKO


m = GEKKO()

maxpop = 8998505 # Maximum population
#iniexpo = 475 # Initally exposed people
#iniinf = 10 # Initially infectd people

# Integration time points:
m.time = np.linspace(0,180,181)

# Constants
L = 271.23 # Influx rate
tcE0 = 3.11e-8 # Transmission between S and E
tcI0 = 0.62e-8 # Transmission between S and I
#tcV0 = m.CV(1e-8) # Transmission between S and V
tcV0 = 1.03e-8 # Transmission between S and V
#tac = m.CV(1e-4) # Transmission adjustment coefficient
tac = 1.01e-4 # Transmission adjustment coefficient


ndr = 3.01e-5  # natural death rate
ddr = 0.01 # disease-induced death rate

ip = 7.   # Incubation period
rr = 1./15. # Recovery rate

rrv = 1. # Removal rate of virus
#srE = m.CV(2.30) # shielding rate by exposed people
srE = 2.30 # shielding rate by exposed people
srI = 0. # shielding rate by infected people


# Population is divided into 4 groups plus the virus concentration (SEIRV):
#S = m.Var(value=89985051000,lb=0,ub=maxpop) # Susceptible
#E = m.Var(value=475,lb=0,ub=maxpop) # Exposed
#I = m.Var(value=10,lb=0,ub=maxpop) # Infected
#R = m.Var(value=10000,lb=0,ub=maxpop) # Recovered
#V = m.Var(value=0,lb=0,ub=maxpop) # Concentration of the virus in the environmental reservoir

pop0 = 8998505

S = m.Var(value=pop0, lb=0, ub=2*pop0) # Susceptible
E = m.Var(value=1000, lb=0, ub=pop0) # Exposed
I = m.Var(value=475, lb=0, ub=pop0) # Infected
R = m.Var(value=10, lb=0, ub=pop0) # Recovered
V = m.Var(value=10000) # Concentration of the virus in the environmental reservoir




# Intermediates:
tcE = m.Intermediate(tcE0/(1.+tac*E))
tcI = m.Intermediate(tcI0/(1.+tac*I))
tcV = m.Intermediate(tcV0/(1.+tac*V))
a = m.Intermediate(1./ip)
X = m.Intermediate(tcE*S*E+tcI*S*I+tcV*S*V)



# Equations:
m.Equation(S.dt()==L-X-ndr*S)
m.Equation(E.dt()==X-(a+ndr)*E)
m.Equation(I.dt()==a*E-(ndr+ddr+rr)*I)
m.Equation(R.dt()==rr*I-ndr*R)
m.Equation(V.dt()==srE*E+srI*I-rrv*V)

m.options.IMODE = 4
m.solve()



# plot results
plt.figure(1)
#plt.plot(m.time,tcE,'b-')
#plt.plot(m.time,tcI,'g--')
#plt.plot(m.time,tcV,'r')


plt.plot(m.time,I,'o-')
plt.plot(m.time,E,'b--')
#plt.plot(m.time,S,'g')
plt.xlabel('Time (days)')
plt.ylabel('Number of people')
plt.legend(['Infected','Exposed', 'Dead'])
plt.show()

