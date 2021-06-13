from pprint import pprint

from astropy import units as u

from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time, TimeDelta
from poliastro.bodies import Earth, Jupiter, Neptune, Saturn, Sun, Uranus
from poliastro.ephem import Ephem
from poliastro.maneuver import Maneuver
from poliastro.plotting import StaticOrbitPlotter, OrbitPlotter2D
from poliastro.twobody import Orbit
from poliastro.util import time_range, norm
from poliastro.threebody import *

solar_system_ephemeris.set("jpl")


# For now, simulate the voyager 2 mission. Have the dates set accordingly for flybys of Jupiter, Saturn. Uranus and Neptune.


launch_date = Time("1977-08-20 14:29", scale="utc").tdb
jupiter_date = Time("1979-07-09 22:29", scale="utc").tdb

# Proper
proper = 1
if proper:
    saturn_date = Time("1981-08-25 03:24", scale="utc").tdb
    uranus_date = Time(
        "1986-01-24", scale="utc"
    ).tdb

    # TODO: Get a precise time, see if that helps at all.
    neptune_date = Time("1989-08-25", scale="utc").tdb
else:
    saturn_date = Time("1983-08-25 03:24", scale="utc").tdb
    uranus_date = Time("1991-01-24", scale="utc").tdb
    neptune_date = Time("2000-08-25", scale="utc").tdb

# 2024
# TODO: Define 2024-ish launch dates and see how they go. First, try to fix the timing problems in the current script

c_3 = 102.4 * u.km ** 2 / u.s ** 2

earth = Ephem.from_body(Earth, time_range(launch_date, end=saturn_date))


plotter = StaticOrbitPlotter()


def assist_to_planet(
    planet,
    target_date,
    current_orbit,
    next_date,
    color="black",
    delta=TimeDelta(0 * 60, format="sec"),
    anom=(0, 180),
    impulses=[]
):
    planet_orbit = Ephem.from_body(planet, time_range(launch_date, end=next_date))
    ss = Orbit.from_ephem(Sun, planet_orbit, target_date)
    man = Maneuver.lambert(current_orbit, ss)
    flyby, target = current_orbit.apply_maneuver(man, intermediate=True)
    for i in impulses:
        pass #TODO
    flyby_end = flyby.propagate(target_date + delta)
    plotter.plot_body_orbit(planet, target_date, label=f"{planet} End of Flyby")
    plotter.plot_trajectory(
        flyby_end.sample(min_anomaly=anom[0] * u.deg, max_anomaly=anom[1] * u.deg),
        label=f"To {planet}",
        color=color,
    )
    return man.impulses, flyby_end, target, man.get_total_cost()


# Basically, we flyby each of the planets successively.

# DONE: Initial launch stuff, and define initial orbit by applying a launch delta-V maneuver to the initial Earth orbit.
ss_e0 = Orbit.from_ephem(Sun, earth, launch_date)
ss_e0_end = ss_e0.propagate_to_anomaly(180.0 * u.deg)
init_orbit = Orbit.from_classical(
    Sun,
    3.76054 * u.AU,
    0.733305558 * u.one,
    4.87003 * u.deg,
    327.8026 * u.deg,
    11.65107 * u.deg,
    348.80709 * u.deg,
)
init_orbit_end = init_orbit.propagate(
    Time("1977-10-20 14:29", scale="utc").tdb
)  # TODO: This bit is sticky.

plotter.plot_body_orbit(Earth, ss_e0_end.epoch, label=f"{Earth} at end of flyby")

jupiter_impulse, jupiter_end, jupiter_target, jupiter_cost = assist_to_planet(
    Jupiter, jupiter_date, init_orbit_end, saturn_date, anom=(0, 150)
)
# plotter.plot_body_orbit(Saturn, saturn_date, label=f"Saturn End of Flyby")
saturn_impulse, saturn_end, saturn_target, saturn_cost = assist_to_planet(
    Saturn, saturn_date, jupiter_end, saturn_date, color="blue", anom=(30, 80)
)
# NOTE: The basic problem seems to be the time.
# To get a curved path, we need a surprising amount of time to get from Jupiter to Saturn.
# I wonder if we can define a single maneuver that takes us on a path from Earth all the way to Saturn.
# I imagine we could specify target orbits rather than simply target points based on ephemeris files,
# but that doesn't work for the next mission.
uranus_impulse, uranus_end, uranus_target, uranus_cost = assist_to_planet(
    Uranus, uranus_date, saturn_end, uranus_date, color="red", anom=(5, 70)
)

neptune_impulse, neptune_end, neptune_target, neptune_cost = assist_to_planet(
    Neptune, neptune_date, uranus_end, neptune_date, color="green", anom=(40, 70)
)

# TODO: Create ephemeris (set of RV Use)
# vectors flyby and target .rv()
def magnitude(obj):
    return (norm(obj.rv()[0]).to(u.AU), norm(obj.rv()[1]))


vectors = [
    magnitude(init_orbit),
    magnitude(jupiter_end),
    magnitude(saturn_end),
    magnitude(uranus_end),
    magnitude(neptune_end),
]

pprint(vectors)
print()
pprint(f"Jupiter: {jupiter_cost}")
print()
pprint(f"Saturn: {saturn_cost}")
print()
pprint(f"Uranus: {uranus_cost}")
print()
pprint(f"Neptune: {neptune_cost}")
print()

print(jupiter_impulse)


# Maneuvers
# Hashtable, with names TCM1, TCM2, etc. Each value is a maneuver with a time specified as an astropy epoch, and a vector of delta-V in
m=Maneuver
tcm={"01":m(()),
     "02":m(()),
     "03":m(()),
     "04":m((epoch,[-563.738,522.016,29.424]*u.mm/u.s)),
     "05":m((epoch,[8723.415,-8807.558,-3908.574]*u.mm/u.s)),
     "06":m((epoch,[-166.348,607.505,500.864]*u.mm/u.s)),
     "07":m(()),
     "08":m((epoch,[-469.649,882.511,648.463]*u.mm/u.s)),
     "09":m((epoch,[-230.733,1205.817,-534.300]*u.mm/u.s)),
     "10":m(()),
     "11":m(()),
     "12":m(()),
     "13":m(()),
     "14":m(()),
     "15":m(()),
     "16":m(()),
     "17":m((epoch,[91.925,-339.117,-10.428]*u.mm/u.s)),
     "18":m((epoch,[-335.093,-850.755,-161.348]*u.mm/u.s)),
     "19":m(()), # Cancelled
     "20":m((epoch,[-477.544,-8.441,-12.920])),
     }
