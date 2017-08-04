import numpy as np
from gpkit import Model, Variable, Vectorize

class Aircraft(Model):
	def setup(self):
		self.wing = Wing()
		self.prop = Propulsion()
		self.components = [self.wing]

		self.mass = Variable('m','kg','mass')
		return self.components, [
			self.mass >= sum(c['m'] for c in self.components)
		]

	def dynamic(self, state):
		return AircraftP(self, state)

class AircraftP(Model):
	def setup(self, aircraft, state):
		self.aircraft = aircraft
		self.wing_aero = aircraft.wing.dynamic(state)
		self.prop_perf = aircraft.prop.dynamic(state)

		self.perf_models = [self.wing_aero,self.prop_perf]

		return self.perf_models, [
			aircraft.mass*state['g'] <= 0.5*state['rho']*(state['V']**2)*self.wing_aero["C_L"]*aircraft.wing['S'],
			self.prop_perf['T'] >= 0.5*state['rho']*(state['V']**2)*self.wing_aero["C_D"]*aircraft.wing['S']
		]

class Wing(Model):
	def setup(self):
		m = Variable('m','kg', "mass")
		S = Variable("S", 190, "m^2", "surface area")
		rho = Variable("\\rho", 1, "kg/m^2", "areal density")
		A = Variable("A", 27, "-", "aspect ratio")
		c = Variable("c", "ft", "mean chord")

		return [m >= S*rho,
		c == (S/A)**0.5]

	def dynamic(self, state):
		return WingAero(self, state)

class WingAero(Model):
    "Wing aerodynamics"
    def setup(self, wing, state):
        CD = Variable("C_D", "-", "drag coefficient")
        CL = Variable("C_L", "-", "lift coefficient")
        e = Variable("e", 0.9, "-", "Oswald efficiency")
        Re = Variable("Re", "-", "Reynold's number")
        D = Variable("D", "lbf", "drag force")

        return [
            CD >= (0.074/Re**0.2 + CL**2/np.pi/wing["A"]/e),
            Re == state["rho"]*state["V"]*wing["c"]/state["\\mu"],
            D >= 0.5*state["rho"]*state["V"]**2*CD*wing["S"],
            ]

class Propulsion(Model):
	def setup(self):
		m = Variable('m','kg','mass')
	def dynamic(self,state):
		return PropPerf(self,state)

class PropPerf(Model):
	def setup(self,prop,state):
		T = Variable('T',1e4,'N','thrust')


class FlightState(Model):
	def setup(self):
		Variable("V", "knots", "true airspeed")
		Variable("\\mu", 1.628e-5, "N*s/m^2", "dynamic viscosity")
		Variable("rho", 0.74, "kg/m^3", "air density")
		Variable('g', 9.8, 'm/s/s', 'acceleration due to gravity')

class FlightSegment(Model):
	def setup(self, aircraft):
		self.flightstate = FlightState()
		self.aircraftp = aircraft.dynamic(self.flightstate)
		return self.flightstate, self.aircraftp

class Mission(Model):
	def setup(self, aircraft):
		with Vectorize(4):
			self.fs = FlightSegment(aircraft)
		self.Vmax = self.fs['V'][3]
		return self.fs

ac = Aircraft()
mission = Mission(ac)
m = Model(1/mission.Vmax,[mission,ac])
sol = m.solve()
print sol.table()