"""Modular aircraft concept"""
import numpy as np
from gpkit import Model, Variable, Vectorize, VectorVariable

class Aircraft(Model):
    "The vehicle model"
    def setup(self):
        self.fuselage = Fuselage()
        self.observable = Observable()

        self.wing = Wing()
        self.flight_deck = FlightDeck()
        self.seat_fl = Seat()
        self.seat_fr = Seat()
        self.seat_rl = Seat()
        self.vtail = Vtail()

        self.turbofan = Turbofan()
        self.gear_pod_left = GearPod()
        self.gear_pod_right = GearPod()

        self.rear_gear_left = RearGear()
        self.rear_gear_right = RearGear()

        self.front_gear = FrontGear()

        # self.pilot = Pilot()
        self.fuel_tank = FuelTank()
        self.fuel = Fuel()
        # self.tail = Tail()
        # self.avionics = Avionics()
        # self.prop_tower = PropTower()
        self.components = [self.hull, self.turbojet, self.drone, self.wing, self.pontoon_left,self.pontoon_right]

        W = Variable("W", "N", "weight")
        self.weight = W
        n_0 = Variable("n_0",0.3*0.7,"-")
        wettedAreaRatio = Variable("wettedAreaRatio","-")
        C_D0 = Variable("C_D0","-")
        return self.components, [
            W >= sum(c["W"] for c in self.components),
            self.wing["S_wet"] >= self.wing["S"]*(1.977),
            wettedAreaRatio >= (self.wing["S_wet"] + self.hull["S_wet"])/self.wing["S"],
            C_D0 >= 0.0065*wettedAreaRatio
        ]

    def dynamic(self, state):
        "This component's performance model for a given state."
        return AircraftP(self, state)

class AircraftP(Model):
    def setup(self, aircraft, state):
        self.aircraft = aircraft
        self.wing_aero = aircraft.wing.dynamic(state)
        # self.engine_p = aircraft.engine.dynamic(state)
        # self.propeller_p = aircraft.propeller.dynamic(state)
        self.turbojet_p = aircraft.turbojet.dynamic(state)
        self.perf_models = [self.wing_aero,self.turbojet_p]
        Wfuel = Variable("W_{fuel}", "N", "fuel weight")
        Wburn = Variable("W_{burn}", "N", "segment fuel burn")
        self.D = Variable("D","N")
        self.Range = Variable("range","m")
        self.z_bre = Variable("z_bre","-")
        self.LD = Variable("LD","-")
        return self.perf_models, [
            aircraft.weight + Wfuel <= (0.5*state["\\rho"]*state["V"]**2
                                        * self.wing_aero["C_L"]
                                        * aircraft.wing["S"]),
            self.D >= sum(p["D"] for p in self.perf_models) + 0.5*state["\\rho"]*state["V"]**2*aircraft["C_D0"]*aircraft.wing["S"],
            Wburn/(state['g']*self.turbojet_p['mdot']) == self.Range/state['V'],
            Wfuel <= state['g']*self.aircraft.fuel_tank['V']*self.aircraft.fuel['rho'], 
            self.turbojet_p["T"] >= self.D,
            self.z_bre >= (state["g"]*self.Range*self.D)/(aircraft.fuel['h']*aircraft["n_0"]*aircraft.weight),
            Wburn/aircraft.weight >= self.z_bre + (self.z_bre**2)/2 + (self.z_bre**3)/6 + (self.z_bre**4)/24,
        ]

class FlightState(Model):
    "Context for evaluating flight physics"
    def setup(self):
        Variable("V","m/s", "true airspeed")
        Variable("\\mu", 1.789e-5, "N*s/m^2", "dynamic viscosity")
        Variable("\\rho", 1.225, "kg/m^3", "air density")
        Variable("g", 9.8, "m/s/s", "gravity")

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
            fs = FlightSegment(aircraft)
        Vmin = Variable('Vmin',10,'m/s')
        Wburn = fs.aircraftp["W_{burn}"]
        Wfuel = fs.aircraftp["W_{fuel}"]
        self.takeoff_fuel = Wfuel[0]
        self.range = Variable("range","m")
        self.z_bre = Variable("z_bre","-")
        return fs, [
                    Wfuel[:-1] >= Wfuel[1:] + Wburn[:-1],
                    Wfuel[-1] >= Wburn[-1],
                    fs.aircraftp.Range[0] == 0.01*self.range,
                    fs.flightstate['V'][0] <= Vmin, 
                    fs.aircraftp.Range[1:] == self.range/3
        ]

class Wing(Model):
    "Aircraft wing model"

    def setup(self):
        W = Variable("W", "N", "weight")
        S = Variable("S", "m^2", "surface area")
        S_wet = Variable("S_wet", "m^2", "wetted area")
        rho = Variable("\\rho", 0.425*9.8, "N/m^2", "areal density")
        A = Variable("A", "-", "aspect ratio")
        c = Variable("c", "m", "mean chord")
        b = Variable("b", "m", "span")
        return [W >= S*rho,
                c == (S/A)**0.5,
                b == A*c,
                S == b*c]
    def dynamic(self, state):
        "Returns this component's performance model for a given state."
        return WingAero(self, state)


class WingAero(Model):
    "Wing aerodynamics"
    def setup(self, wing, state):
        C_D = Variable("C_D", "-", "drag coefficient")
        CL = Variable("C_L", "-", "lift coefficient")
        e = Variable("e", "-", "Oswald efficiency")
        Re = Variable("Re", "-", "Reynold's number")
        D = Variable("D", "N", "drag force")
        sigma = Variable("sigma","-","ground effect correction factor")  #p148 YBD
        return [
            e + (1.78*0.045*wing["A"]**0.68) <= 1.14, 
            e >= 0.1,
            wing["A"]<=3*1.38, #Raymer 12.6
            C_D >= (0.074/(Re**0.2) + CL**2/(np.pi*wing["A"]*e)),
            Re == state["\\rho"]*state["V"]*wing["c"]/state["\\mu"],
            D >= 0.5*state["\\rho"]*state["V"]**2*C_D*wing["S"]
        ]


class Fuselage(Model):
    "The thing that carries the fuel, engine, and payload"
    def setup(self):
        m = Variable("m", "kg", "mass")
        S_wet = Variable("S_wet", 3, "m^2")
    def dynamic(self,hull,state):
        return FuselageP(self,state)

class FuselageP(Model):
    def setup(self,fuselage,state):
        m = Variable("m","kg","mass")

class Observable(Model):
    def setup(self):
        m = Variable("m","kg","mass")

class Seat(Model):
    def setup(self):
        m = Variable("m",32,"kg","mass") #Raymer pg 456

class FlightDeck(Model):
    def setup(self):
        m = Variable("m",60,"kg","mass") #Raymer pg 456

class Vtail(Model):
    def setup(self):
        m = Variable("m",)

# class Engine(Model):
#     def setup(self):
#         W = Variable("W",421.4,"N","weight")
#         P_max = Variable("P_max",26099,"W","max engine power output")
#         BSFC_min = Variable("BSFC_min",1.4359644493044238e-07,"kg/W/s","maximum brake specfic fuel consumption")
#     def dynamic(self,state):
#         return EngineP(self,state)

# class EngineP(Model):
#     def setup(self,engine,state):
#         mdot = Variable("mdot","kg/s","fuel mass flow")
#         self.P = Variable("P","W","engine power output")
#         throttle = Variable("throttle","-","throttle position")
#         D = Variable("D",1e-7,"N")
#         return [self.P <= mdot/engine["BSFC_min"],
#                 self.P <= engine["P_max"],
#                 throttle == self.P/engine["P_max"]]

# class Propeller(Model):
#     def setup(self):
#         W = Variable("W",50,"N","weight")
#         diameter = Variable("diameter",1.5,"m","Diameter")
#     def dynamic(self,state):
#         return PropellerP(self,state)

# class PropellerP(Model):
#     def setup(self,propeller,state):
#         # What the hell is a?
#         a = 0.8
#         T = Variable("T","N","thrust")
#         f = Variable("f","1/s","frequency")
#         P_in = Variable("P_in","W","power input")
#         n = Variable("n",0.8,"-","efficiency")
#         J = Variable("J","-","advance ratio")
#         D = Variable("D",1e-7,"N")
#         return [J == state['V']/(f*propeller["diameter"]),
#                 f <= Variable("f",45,"1/s","prop speed limit"),
#                 n <= a*J,
#                 T<= P_in*n/state['V']]

class Turbojet(Model):
    def setup(self):
        # Specs for PBS TJ20
        W = Variable("W",2.1*9.8,"N")
        T_max = Variable("T_max",210,"N")
        TSFC = Variable("TSFC",0.000277778*0.165,"kg/(N*s)")
    def dynamic(self,state):
        return TurbojetP(self,state)

class TurbojetP(Model):
    def setup(self,turbojet,state):
        T = Variable("T","N")
        mdot = Variable("mdot","kg/s")
        D = Variable("D",1e-7,"N")
        return [    T<=turbojet["T_max"],
                    T<=mdot/turbojet["TSFC"]]

class Pilot(Model):
    def setup(self):
        W = Variable("W",735,"N","weight")

class FuelTank(Model):
    def setup(self):
        V = Variable("V","m^3","volume")
        l = Variable("l","m")
        w = Variable("w",0.3,"m")
        h = Variable("h",0.2,"m")
        return [V == l*w*h]

class Fuel(Model):
    def setup(self):
        # Typical values for gasoline
        h = Variable("h",42.448e6,"J/kg","heating value")
        rho = Variable("rho",719.7,"kg/m^3","density")

class Tail(Model):
    def setup(self):
        W = Variable("W",98,"N","weight")

class Avionics(Model):
    def setup(self):
        W = Variable("W",98,"N","weight")

class PropTower(Model):
    def setup(self):
        W = Variable("W",16*9.8,"N","weight")

class Drone(Model):
    def setup(self):
        W = Variable("W",4*9.8,"N","maximum weight of inspire 2")

class Pontoon(Model):
    def setup(self):
        W = Variable("W",2*9.8,"N","weight of bouyancy pontoon")


AC = Aircraft()
MISSION = Mission(AC)
# objective = MISSION.range[0] + MISSION.range[1] + MISSION.range[2] + MISSION.range[3]
M = Model(1/MISSION.range, [MISSION, AC])
M.debug()
SOL = M.solve(verbosity=0)
print SOL.table()
rangeInKm = 0.001*SOL(MISSION.range)
#We want the range to be 500 km, as determined by the crossing from Greenland to Iceland and Iceland to Faroe Islands
#these crossings are around 480 km, which means we need to have an ideal cruise in that wave height of 500 km
#without that, no point building the thing
rangeObjective = 500

print('Range is ' + str(rangeInKm)+' km')
print('Fuel tank length is ' + str(SOL(AC.fuel_tank["V"])/0.6))
# print('At %s percent of range goal'%str(100*rangeInKm/rangeObjective))
# print('MTOW of ' + str(sol(W)))
# print('Cruise LD of '+str(sol(C_L/C_D)))
# W_fuelsol = sol(MISSION[])
# W_zfwsol = sol(W_zfw)
# print('Full fuel fraction of ' +str(W_fuelsol/(W_zfwsol+W_fuelsol)))