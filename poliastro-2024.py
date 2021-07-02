import warnings

warnings.filterwarnings("ignore")

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


def magnitude(obj):
    return (norm(obj.rv()[0]).to(u.AU), norm(obj.rv()[1]))


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


solar_system_ephemeris.set("jpl")
C_3 = 102.4 * u.km ** 2 / u.s ** 2

dates = {
    "launch": Time("2024-01-01", scale="utc").tdb,
    "jupiter": Time("2026-01-01", scale="utc").tdb,
    "uranus": Time("2033-01-01", scale="utc").tdb,
    "neptune": Time("2036-01-01", scale="utc").tdb,
}
for i in range(20):
    for k in dates.keys():
        dates[k] = dates[k] + TimeDelta(12 * u.wk)

    plotter = StaticOrbitPlotter()
    earth = Ephem.from_body(Earth, time_range(dates["launch"], end=dates["jupiter"]))

    ss_e0 = Orbit.from_ephem(Sun, earth, dates["launch"])
    ss_e0_end = ss_e0.propagate_to_anomaly(180.0 * u.deg)
    init = Orbit.synchronous(Earth).change_attractor(Sun, force=True)
    plotter.plot_body_orbit(Earth, ss_e0_end.epoch, label=f"{Earth} at end of flyby")

    ends = {"init": init.propagate(Time("2024-12-20 14:29", scale="utc").tdb)}
    costs = {}

    ends["jupiter"], costs["jupiter"] = assist_to_planet(
        Jupiter, dates["jupiter"], ends["init"], "red"
    )
    ends["uranus"], costs["uranus"] = assist_to_planet(
        Uranus, dates["uranus"], ends["jupiter"], "green"
    )
    ends["neptune"], costs["neptune"] = assist_to_planet(
        Neptune, dates["neptune"], ends["uranus"], "blue"
    )

for i in ("Jupiter", "Uranus", "Neptune"):
    pprint(f"{i}: {costs[i.lower()]}")
pprint([magnitude(i) for i in ends.values()])
