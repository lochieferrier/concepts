import numpy as np
from gpkit import Model, Variable, Vectorize

##Flight segments
## 0 - Nothing (climb to record altitude?)
## 1 - 60 second cruise at level condition
## 2 - The acceleration period to be later simulated
## 3 - 30 second supercruise

class Aircraft(Model):
	def setup(self):
		self.wing = Wing()
		self.prop = Propulsion()
		self.battery = Battery()
		self.components = [self.wing,self.prop,self.battery]

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
		self.acceleration = Variable('a','m/s/s','acceleration or thrust excess')

		return self.perf_models, [
			aircraft.mass*state['g'] <= 0.5*state['rho']*(state['V']**2)*self.wing_aero["C_L"]*aircraft.wing['S'],
			self.prop_perf['T'] >= self.wing_aero['D']
		]

class Wing(Model):
	def setup(self):
		m = Variable('m','kg', "mass")
		S = Variable("S","m^2", "surface area")
		rho = Variable("\\rho", 1, "kg/m^2", "areal density")
		A = Variable("A","-", "aspect ratio")
		c = Variable("c", "ft", "mean chord")
		C_Lmax = Variable('C_Lmax',1.2,'-','maximum design C_L')
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
            CL <= wing['C_Lmax']
        ]
class Battery(Model):
	def setup(self):
		m = Variable('m','kg','mass')
		E = Variable('E','joule','battery capacity')
		Estar = Variable('Estar',300,'watt_hour/kg','battery pack specific energy')
		constraints = [m >= E/Estar]
		return constraints

	def dynamic(self, state):
		return BatteryPerf(self,state)

class BatteryPerf(Model):
	def setup(self,battery,state):
		V = Variable('V','V','battery voltage')

class Propulsion(Model):
	def setup(self):
		m = Variable('m','kg','mass')
		Pstar = Variable('Pstar',1e3,'W/kg','relationship between shaft power and propulsion mass')
		Pmax = Variable('Pmax','W','Maximum sustainable shaft power passable by motor')
		constraints = [m >= Pmax/Pstar]
		return constraints

	def dynamic(self,state):
		return PropPerf(self,state)

class PropPerf(Model):
	def setup(self,prop,state):
		T = Variable('T','N','thrust')
		P = Variable('P','W','power')
		# 'Efficiency' stat derived from the Olympus 593
		# 39 g / s for a kN, and a gram of Jet A-1 is 
		Estar_A1 = Variable('Estar_A1',42.8e6,'J/kg')
		SFC_comparative = Variable('SFC_comparative',33.8,'g/(s*kN)')
		n = Variable('n','N/W','efficiency')
		return [T <= P*n,
				P <= prop['Pmax'],
				n == 1/(SFC_comparative*Estar_A1)]

class FlightState(Model):
	def setup(self):
		V = Variable("V", "m/s", "true airspeed")
		h = Variable('h','ft','altitude')
		rho = Variable('rho','kg/m^3','air density')
		# VERY BAD, JUST TAKING THE MIDDLE OF THE RANGE FOR MU
		Variable("\\mu", 8.45e-5, "N*s/m^2", "dynamic viscosity")
		Variable('g', 9.8, 'm/s/s', 'acceleration due to gravity')
		a = Variable('a',295.1,'m/s','speed of sound')
		M = Variable('M','-','mach number the aircraft is flying at')
		constraints = [M == V/a]
		return constraints
		# # self.atm = Atmosphere(h)
		# return self.atm

# class Atmosphere(Model):
#     """
#     Model to capture density changes with altitude
#     """
#     def __init__(self, h,N=1, **kwargs):

#         # h = Variable("h", alt, "ft", "Altitude")
#         p_sl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
#         L_atm = Variable("L_{atm}", 0.0065, "K/m", "Temperature lapse rate")
#         T_sl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")

#         T_atm = Variable("T_{atm}","K", "Air temperature")
#         mu_atm = Variable("\\mu", "N*s/m^2", "Dynamic viscosity")
#         mu_sl = Variable("\\mu_{sl}", 1.789*10**-5, "N*s/m^2",
#                          "Dynamic viscosity at sea level")
#         R_spec = Variable("R_{spec}", 287.058, "J/kg/K",
#                           "Specific gas constant of air")
#         h_ref = Variable("h_{ref}", 15000, "ft", "Reference altitude")
#         rho = Variable("\\rho", "kg/m^3", "Air density")

#         # Atmospheric variation with altitude (valid from 0-7km of altitude)
#         constraints = [T_sl >= T_atm + L_atm*h ,
#         			   rho == p_sl*T_atm**(5.257-1)/R_spec/(T_sl**5.257),
#                        (mu_atm/mu_sl)**0.1 == 0.991*(h/h_ref)**(-0.00529)
#                       ]

#         Model.__init__(self, None, constraints, **kwargs)

class FlightSegment(Model):
	def setup(self, aircraft):
		self.flightstate = FlightState()
		self.aircraftp = aircraft.dynamic(self.flightstate)
		self.t = Variable('t','s','duration of flight segment')
		self.E_used = Variable('E_used','joule','energy used from battery in flight segment')

		return self.flightstate, self.aircraftp, [self.E_used >= self.aircraftp.prop_perf['P']*self.t]

class Mission(Model):
	def setup(self, aircraft):
		subsonicFlatCruiseTime = Variable('t_sub_flat',60,'second')
		supersonicCruiseTime = Variable('t_super_flat',30,'second')
		with Vectorize(3):
			self.fs = FlightSegment(aircraft)
		E_used = self.fs['E_used']
		self.Vmax = self.fs['V'][1]
		self.Vmin = Variable('Vmin',100,'m/s','stall speed for a reasonable landing')
		self.MTOW = aircraft.mass
		return self.fs,[aircraft.battery['E'] >= sum(E for E in E_used),
		self.fs.t[0] >= subsonicFlatCruiseTime,
		self.fs.t[1] >= supersonicCruiseTime,
		self.fs.t[2] >= Variable('landtime',10,'s'),
		self.fs['rho'][0] == self.fs['rho'][1],
		self.fs['V'][2] <= self.Vmin,
		self.fs['rho'][0] >= Variable('rho_llim',0.07258,'kg/m^3'),
		self.fs['rho'][0] <= Variable('rho_ulim',0.2755,'kg/m^3'),
		self.fs['rho'][2] == Variable('rho_sea',1.225,'kg/m^3'),
		aircraft.mass <= Variable('masslim',25,'kg'),
		self.fs['M'][1] >= Variable('supersonic_req',1,'-')
		# self.fs[0].aircraftp.acceleration
	]

ac = Aircraft()
mission = Mission(ac)
m = Model(1/mission.Vmax,[mission,ac])
m.debug()
sol = m.solve()

print sol.table()