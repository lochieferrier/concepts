import numpy as np
from gpkit import Model, Variable, VectorVariable

#Determine the best design for a max range rocket over 10 timesteps
N = 10
m = VectorVariable(N,'m','kg','total mass')
V = VectorVariable(N,'V','m/s','velocity')
dV = VectorVariable(N,'dV',[2]*N,'m/s','difference in velocity over a segment')
x = VectorVariable(N,'x','m','total distance flown')
S = Variable(2,'S','m^2','Drag reference area')

constraints = [V[N-1] >= sum(dV for dV in V[:N-1])]
# print constraints
objective = 1/V[N-1]
print objective,constraints
m = Model(constraints,objective)
