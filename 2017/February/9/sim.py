"""Modular aircraft concept"""
import numpy as np
from gpkit import Model, Variable, Vectorize


class Aircraft(Model):
    "The vehicle model"
    def setup(self):
        self.fuse = Fuselage()
        self.wing = Wing()
        self.drone = Drone()
        self.fuel = Fuel()
        self.components = [self.fuse, self.wing, self.drone]
        n_0 = Variable("n_0",0.3*0.7,"-")
        W = Variable("W", "N", "weight")
        self.weight = W

        return self.components, [
            W >= sum(c["W"] for c in self.components)
        ]

    def dynamic(self, state):
        "This component's performance model for a given state."
        return AircraftP(self, state)


class AircraftP(Model):
    "Aircraft flight physics: weight <= lift, fuel burn"
    def setup(self, aircraft, state):
        self.aircraft = aircraft
        self.wing_aero = aircraft.wing.dynamic(state)
        self.perf_models = [self.wing_aero]
        Wfuel = Variable("W_{fuel}", "N", "fuel weight")
        Wburn = Variable("W_{burn}", "N", "segment fuel burn")
        self.Range = Variable("range","m")
        self.z_bre = Variable("z_bre","-")
        self.D = Variable("D","N")

        return self.perf_models, [
            self.D >= self.wing_aero["D"],
            aircraft.weight + Wfuel <= (0.5*state["\\rho"]*state["V"]**2
                                        * self.wing_aero["C_L"]
                                        * aircraft.wing["S"]),
            Wburn/(state['g']*state['m_dot']) == self.Range/state['V'],
            self.z_bre >= (state["g"]*self.Range*self.D)/(aircraft.fuel['h']*aircraft["n_0"]*aircraft.weight),
            Wburn/aircraft.weight >= self.z_bre + (self.z_bre**2)/2 + (self.z_bre**3)/6 + (self.z_bre**4)/24,
    ]


class FlightState(Model):
    "Context for evaluating flight physics"
    def setup(self):
        Variable("V", 20, "m/s", "true airspeed")
        Variable("m_dot","kg/s","fuel flow")
        Variable("\\mu", 1.73e-5, "N*s/m^2", "dynamic viscosity")
        Variable("\\rho", 1.225, "kg/m^3", "air density")
        Variable("g",9.8,"m/s/s")

class FlightSegment(Model):
    "Combines a context (flight state) and a component (the aircraft)"
    def setup(self, aircraft):
        self.flightstate = FlightState()
        self.aircraftp = aircraft.dynamic(self.flightstate)
        return self.flightstate, self.aircraftp


class Mission(Model):
    "A sequence of flight segments"
    def setup(self, aircraft):
        with Vectorize(4):  # four flight segments
            self.fs = FlightSegment(aircraft)

        Wburn = self.fs.aircraftp["W_{burn}"]
        Wfuel = self.fs.aircraftp["W_{fuel}"]
        self.takeoff_fuel = Wfuel[0]

        return self.fs, [Wfuel[:-1] >= Wfuel[1:] + Wburn[:-1],
                         Wfuel[-1] >= Wburn[-1]]


class Wing(Model):
    "Aircraft wing model"
    def dynamic(self, state):
        "Returns this component's performance model for a given state."
        return WingAero(self, state)

    def setup(self):
        W = Variable("W", "lbf", "weight")
        S = Variable("S", 190, "ft^2", "surface area")
        rho = Variable("\\rho", 1, "lbf/ft^2", "areal density")
        A = Variable("A", 27, "-", "aspect ratio")
        c = Variable("c", "ft", "mean chord")

        return [W >= S*rho,
                c == (S/A)**0.5]


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
            Re == state["\\rho"]*state["V"]*wing["c"]/state["\\mu"],
            D >= 0.5*state["\\rho"]*state["V"]**2*CD*wing["S"],
            ]

class Drone(Model):
    def setup(self):
        W = Variable("W",4*9.8,"N")

class Fuselage(Model):
    "The thing that carries the fuel, engine, and payload"
    def setup(self):
        # fuselage needs an external dynamic drag model,
        # left as an exercise for the reader
        # V = Variable("V", 16, "gal", "volume")
        # d = Variable("d", 12, "in", "diameter")
        # S = Variable("S", "ft^2", "wetted area")
        # cd = Variable("c_d", .0047, "-", "drag coefficient")
        # CDA = Variable("CDA", "ft^2", "drag area")
        Variable("W", 100, "lbf", "weight")

class Engine(Model):
    def setup(self):
        # Specs based on PBS TJ20
        W = Variable("W",2.1*9.8,"N")
        T_max = Variable("T_max",210,"N")
        TSFC = Variable("TSFC",0.000277778*0.165,"kg/(N*s)")
    def dynamic(self,state):
        return EngineP(self,state)

class EngineP(Model):
    def setup(self,engine,state):
        T = Variable("T","N")
        return [    T<=T_max,
                    T<=engine["TSFC"]*state["m_dot"]]

class Fuel(Model):
    def static(self):
        h = Variable("h",43.15e6,"J/kg")
        rho = Variable("rho",804,"kg/m^3")

AC = Aircraft()
MISSION = Mission(AC)
M = Model(MISSION.takeoff_fuel, [MISSION, AC])
sol = M.solve(verbosity=0)
print sol.table()