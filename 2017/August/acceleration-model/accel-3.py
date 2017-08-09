"Very simple problem: minimize x while keeping x greater than 1."
from gpkit import Variable, Model,Vectorize

# Decision variable
dV = Variable('dV','m/s')
N=2
with Vectorize(N):
	a = Variable('a',2,'m/s/s')
	t = Variable('t',1,'s')
# Constraint
constraints = [dV == a*t]

# Objective (to minimize)
objective = dV

# Formulate the Model
m = Model(objective, constraints)

# Solve the Model
sol = m.solve(verbosity=0)

print sol.table()