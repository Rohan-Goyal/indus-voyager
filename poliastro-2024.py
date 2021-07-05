import warnings


from pprint import pprint

from astropy import units as u
from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time, TimeDelta
from poliastro.bodies import Venus, Earth, Jupiter, Neptune, Saturn, Sun, Uranus
from poliastro.ephem import Ephem
from poliastro.maneuver import Maneuver
from poliastro.plotting import StaticOrbitPlotter
from poliastro.twobody import Orbit
from poliastro.util import time_range, norm
from numpy import arccos, pi


def magnitude(obj):
    return (norm(obj.rv()[0]).to(u.AU), norm(obj.rv()[1]))


def position(body, epoch):
    return (
        Ephem.from_body(body, epoch).rv()[0]
        if type(epoch) == Time
        else Ephem.from_body(body, Time(epoch)).rv()[0]
    )


def phasing_angle(body1, body2, epoch1, epoch2):
    # Angle of orbit deflections between body1 and body2
    pos1 = position(body1, epoch1)
    pos2 = position(body2, epoch2)
    x = (norm(pos1) ** 2 + norm(pos2 - pos1) ** 2 - norm(pos2) ** 2) / (
        2 * norm(pos1) * norm(pos2 - pos1)
    )
    return arccos(x).to(u.deg)


def angle(body1, body2, epoch1, epoch2):
    # Angle from sun between body1 and body2
    pos1 = position(body1, epoch1)
    pos2 = position(body2, epoch2)
    x = (norm(pos1) ** 2 + norm(pos2) ** 2 - norm(pos2 - pos1) ** 2) / (
        2 * norm(pos1) * norm(pos2)
    )
    return arccos(x).to(u.deg)


def jup_phase(launch):
    return angle(Jupiter, Uranus, launch + 2 * u.yr, launch + 9.5 * u.yr)


def loop():  # Iterate over possible launch times, view angle between jupiter and uranus
    launch = Time("2024-01-01")
    for i in range(0, 20):
        l = launch + TimeDelta(i * u.yr)
        print(l, jup_phase(l))


# NOTE: Launch between 2028 and 2033 produces values <55


def assist_to_planet(planet, arrival, current, color):
    planet_orbit = Ephem.from_body(
        planet, time_range(dates["launch"], end=arrival + TimeDelta(1 * u.yr))
    )
    ss = Orbit.from_ephem(Sun, planet_orbit, arrival)
    man = Maneuver.lambert(current, ss)
    (flyby, _) = current.apply_maneuver(man, intermediate=True)
    flyby_end = flyby.propagate(arrival + TimeDelta(0 * u.wk))
    plotter.plot_body_orbit(planet, arrival, label=f"{planet} Time of Flyby")
    plotter.plot_trajectory(
        flyby_end.sample(min_anomaly=flyby.nu, max_anomaly=flyby_end.nu),
        label=f"To {planet}",
        color=color,
    )
    return (
        flyby_end,
        man.get_total_cost(),
    )


warnings.filterwarnings("ignore")
solar_system_ephemeris.set("jpl")
C_3 = 102.4 * u.km ** 2 / u.s ** 2

launchdate = Time("2029-04-01", scale="utc").tdb
dates = {
    "launch": launchdate,
    "jupiter": launchdate + TimeDelta(2 * u.yr),
    "uranus": launchdate + TimeDelta(9.5 * u.yr),
    "neptune": launchdate + TimeDelta(12 * u.yr),
}

plotter = StaticOrbitPlotter()
earth = Ephem.from_body(Earth, time_range(dates["launch"], end=dates["jupiter"]))

ss_e0 = Orbit.from_ephem(Sun, earth, dates["launch"])
ss_e0_end = ss_e0.propagate_to_anomaly(180.0 * u.deg)
init = Orbit.synchronous(Earth, epoch=dates["launch"]).change_attractor(Sun, force=True)
plotter.plot_body_orbit(Earth, ss_e0_end.epoch, label=f"{Earth} at end of flyby")

ends = {"init": init.propagate(Time("2026-12-20 14:29", scale="utc").tdb)}
costs = {}

ends["jupiter"], costs["jupiter"] = assist_to_planet(
    Jupiter, dates["jupiter"], ends["init"], "red"
)
ends["uranus"], costs["uranus"] = assist_to_planet(
    Uranus, dates["uranus"], ends["jupiter"], "green"
)
# ends["neptune"], costs["neptune"] = assist_to_planet(
#     Neptune, dates["neptune"], ends["uranus"], "blue"
# )

pprint([f"{i.capitalize()}: {v} \n" for i, v in costs.items()])
pprint([f"{i.capitalize()}:{magnitude(v)}" for i, v in ends.items()])
