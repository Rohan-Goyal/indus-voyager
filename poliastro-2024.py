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

plotter = StaticOrbitPlotter()

def magnitude(obj):
    return (norm(obj.rv()[0]).to(u.AU), norm(obj.rv()[1]))


def assist_to_planet(planet, arrival, current, color, anom):
    planet_orbit = Ephem.from_body(
        planet, time_range(launch_date, end=arrival + TimeDelta(1 * u.yr))
    )
    ss = Orbit.from_ephem(Sun, planet_orbit, arrival)
    man = Maneuver.lambert(current, ss)
    (flyby, _) = current.apply_maneuver(man, intermediate=True)
    flyby_end = flyby.propagate(arrival + TimeDelta(0 * u.wk))
    plotter.plot_body_orbit(planet, arrival, label=f"{planet} Time of Flyby")
    plotter.plot_trajectory(
        flyby_end.sample(min_anomaly=anom[0] * u.deg, max_anomaly=anom[1] * u.deg),
        label=f"To {planet}",
        color=color,
    )
    return (
        flyby_end,
        man.get_total_cost(),
    )


launch_date = Time("2024-08-20 14:29", scale="utc").tdb
jupiter_date = Time("2027-07-09 22:29", scale="utc").tdb
uranus_date = Time("2033-07-24", scale="utc").tdb
neptune_date = Time("2036-12-25", scale="utc").tdb

earth = Ephem.from_body(Earth, time_range(launch_date, end=jupiter_date))
C_3 = 102.4 * u.km ** 2 / u.s ** 2
# dv = C_3 ** 0.5 * v_e0 / norm(v_e0)
# man = Maneuver.impulse(dv)

# DONE: Initial launch stuff, and define initial orbit by applying a launch delta-V maneuver to the initial Earth orbit.
ss_e0 = Orbit.from_ephem(Sun, earth, launch_date)
ss_e0_end = ss_e0.propagate_to_anomaly(180.0 * u.deg)
init_orbit = Orbit.synchronous(Earth).change_attractor(Sun, force=True)

solar_system_ephemeris.set("jpl")
# init_orbit=ss_e0.apply_maneuver(man)

# https://astronomy.stackexchange.com/questions/4823/finding-the-radius-of-an-eccentric-orbit-at-any-point

init_orbit_end = init_orbit.propagate(
    Time("2024-12-20 14:29", scale="utc").tdb
)
plotter.plot_body_orbit(Earth, ss_e0_end.epoch, label=f"{Earth} at end of flyby")
plotter.plot_trajectory(init_orbit.sample(min_anomaly=0*u.deg,max_anomaly=180*u.deg))

jupiter_end, jupiter_cost = assist_to_planet(
    Jupiter, jupiter_date, init_orbit_end, "green", (-100,250)
)
uranus_end, uranus_cost = assist_to_planet(
    Uranus, uranus_date, jupiter_end, "red", (-60, 156)
)
# neptune_end, neptune_cost = assist_to_planet(
#     Neptune, neptune_date, uranus_end, "green", (40, 68)
# )

vectors = [
    magnitude(init_orbit),
    magnitude(jupiter_end),
    magnitude(uranus_end),
    # magnitude(neptune_end),
]

pprint(vectors)
pprint(f"Jupiter: {jupiter_cost}")
pprint(f"Uranus: {uranus_cost}")
# pprint(f"Neptune: {neptune_cost}")
