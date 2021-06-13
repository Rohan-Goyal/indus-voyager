from pprint import pprint

from astropy import units as u
from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time, TimeDelta
from poliastro.bodies import Earth, Jupiter, Neptune, Saturn, Sun, Uranus
from poliastro.ephem import Ephem
from poliastro.maneuver import Maneuver
from poliastro.plotting import StaticOrbitPlotter
from poliastro.twobody import Orbit
from poliastro.util import time_range

solar_system_ephemeris.set("jpl")

# TODO
launch_date = Time("2024-08-20 14:29", scale="utc").tdb
jupiter_date = Time("2026-07-09 22:29", scale="utc").tdb
uranus_date = Time("2033-01-24", scale="utc").tdb
neptune_date = Time("2036-08-25", scale="utc").tdb


c_3 = 102.4 * u.km ** 2 / u.s ** 2

earth = Ephem.from_body(Earth, time_range(launch_date, end=jupiter_date))


def assist_to_planet(
    planet,
    target_date,
    current_orbit,
    next_date,
    color="C1",
    delta=TimeDelta(0 * 60, format="sec"),
    anom=(0, 180),
):
    planet_orbit = Ephem.from_body(planet, time_range(launch_date, end=next_date))
    ss = Orbit.from_ephem(Sun, planet_orbit, target_date)
    man = Maneuver.lambert(current_orbit, ss)
    flyby, target = current_orbit.apply_maneuver(man, intermediate=True)
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
init_orbit=Orbit.geostationary(Earth, period=24*u.hour)
init_orbit_end = init_orbit.propagate(
    Time("2024-12-20 14:29", scale="utc").tdb
)  # TODO: This bit is sticky.
plotter = StaticOrbitPlotter()
plotter.plot_body_orbit(Earth, ss_e0_end.epoch, label=f"{Earth} at end of flyby")

jupiter_impulse, jupiter_end, jupiter_target, jupiter_cost = assist_to_planet(
    Jupiter, jupiter_date, init_orbit_end, uranus_date, anom=(0, 150)
)
pprint(f"Jupiter: {jupiter_cost}")
print()

uranus_impulse, uranus_end, uranus_target, uranus_cost = assist_to_planet(
    Uranus, uranus_date, jupiter_end, neptune_date, color="C3", anom=(20, 95)
)
pprint(f"Uranus: {uranus_cost}")
print()

neptune_impulse, neptune_end, neptune_target, neptune_cost = assist_to_planet(
    Neptune, neptune_date, uranus_end, neptune_date, color="C4", anom=(80, 120)
)
pprint(f"Neptune: {neptune_cost}")
print()
