from pprint import pprint
from astropy import units as u

from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time, TimeDelta
from poliastro.bodies import Earth, Jupiter, Neptune, Saturn, Sun, Uranus
from poliastro.ephem import Ephem
from poliastro.maneuver import Maneuver
from poliastro.plotting import StaticOrbitPlotter
from poliastro.twobody import Orbit
from poliastro.util import time_range, norm


plotter = StaticOrbitPlotter()


def magnitude(obj):
    return (norm(obj.rv()[0]).to(u.AU), norm(obj.rv()[1]))


def gen_ephem(orbit, start):
    # Output RV vectors and epochs.
    # One implementation is simply propagate to time + timedelta*i, output epoch and rv to a file.
    # Another is more complex, but probably cheaper. Basically, start by taking the sample from a1 to a2.
    # TODO: Picking A1 and A2. We want readings every few minutes, over a 24-hour period (I guess).
    # We want every 10 minutes (say, meaning a total of 144 samples over 24 hours).
    # I'm stupid. It's quite simple: Find anom at 00:00 on date, and find anom at 23:59 on date. Then sample
    init = orbit.propagate(start)
    final = orbit.propagate(start + TimeDelta(1 * u.d))
    data = orbit.sample(values=144, init.nu, final.nu)  # Every 10 minutes.
    return data, data.differentials["s"]  # Position in km, vel in km/s


def gen_ephem_backup(orbit, start, step=10 * 60 * u.s):
    return [
        (start + i * step, orbit.propagate(start + i * step).rv()) for i in range(144)
    ]


def tcm():

    t = Time
    m = Maneuver.impulse

    # Maneuvers
    # Hashtable, with names TCM1, TCM2, etc. Each value is a maneuver with a time specified as an astropy epoch, and a vector of delta-V in
    return {
        # "01":m(()),
        # "02":m(()),
        # "03":m(()),
        "04": (
            m([563.738, 522.016, 29.424] * u.mm / u.s),
            t("1979-06-27 10:09:09", scale="utc").tdb,
        ),
        "05": (
            m([8723.415, -8807.558, -3908.574] * u.mm / u.s),
            t("1979-07-10 00:39:39", scale="utc").tdb,
        ),
        "06": (
            m([-166.348, 607.505, 500.864] * u.mm / u.s),
            t("1979-07-23 16:09:50", scale="utc").tdb,
        ),
        # "07":m(()),
        "08": (
            m([-469.649, 882.511, 648.463] * u.mm / u.s),
            t("1981-07-19 11:16:25", scale="utc").tdb,
        ),
        "09": (
            m([-230.733, 1205.817, -534.300] * u.mm / u.s),
            t("1981-08-18 21:26:16", scale="utc").tdb,
        ),
        # "10":m(()),
        # "11":m(()),
        # "12":m(()),
        # "13":m(()),
        # "14":m(()),
        # "15":m(()),
        # "16":m(()),
        "17": (
            m([91.925, -339.117, -10.428] * u.mm / u.s),
            t("1989-04-20 16:19:46", scale="utc").tdb,
        ),
        "18": (
            m([-335.093, -850.755, -161.348] * u.mm / u.s),
            t("1989-08-01 12:55:18", scale="utc").tdb,
        ),
        # "19":m(()), # Cancelled
        "20": (
            m([-477.544, -8.441, -12.920] * u.mm / u.s),
            t("1989-08-21 12:48:06", scale="utc").tdb,
        ),
    }


def impulses(date, current, mans=tcm()):
    # Get all impulses between end of current orbit and target dat
    return [i for i in mans.values() if current.epoch < i[1] < date]


def assist_to_planet(planet, arrival, current, color, anom):
    planet_orbit = Ephem.from_body(
        planet, time_range(launch_date, end=arrival + TimeDelta(1 * u.yr))
    )
    ss = Orbit.from_ephem(Sun, planet_orbit, arrival)
    man = Maneuver.lambert(current, ss)
    flyby, target = current.apply_maneuver(man, intermediate=True)
    for i in impulses(arrival, current):
        flyby.propagate(i[1])
        (flyby,) = flyby.apply_maneuver(i[0], intermediate=True)
    flyby_end = flyby.propagate(arrival + TimeDelta(0 * u.wk))
    plotter.plot_body_orbit(planet, arrival, label=f"{planet} End of Flyby")
    plotter.plot_trajectory(
        flyby_end.sample(min_anomaly=anom[0] * u.deg, max_anomaly=anom[1] * u.deg),
        label=f"To {planet}",
        color=color,
    )
    return (
        flyby_end,
        man.get_total_cost(),
        gen_ephem(flyby, Time(arrival.iso.split()[0] + " 00:00")),
    )


# Removed returns of man.impulses and target


# For now, simulate the voyager 2 mission. Have the dates set accordingly for flybys of Jupiter, Saturn. Uranus and Neptune.
solar_system_ephemeris.set("jpl")

launch_date = Time("1977-08-20 14:29", scale="utc").tdb
jupiter_date = Time("1979-07-09 22:29", scale="utc").tdb
saturn_date = Time("1981-08-25 03:24", scale="utc").tdb
uranus_date = Time("1986-01-24 17:59", scale="utc").tdb
neptune_date = Time("1989-08-25 03:56", scale="utc").tdb

earth = Ephem.from_body(Earth, time_range(launch_date, end=saturn_date))


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
    # epoch=launch_date
)
init_orbit_end = init_orbit.propagate(
    Time("1977-10-20 14:29", scale="utc").tdb
)  # TODO: This bit is sticky.

plotter.plot_body_orbit(Earth, ss_e0_end.epoch, label=f"{Earth} at end of flyby")


# NOTE: The basic problem seems to be the time.
# To get a curved path, we need a surprising amount of time to get from Jupiter to Saturn.
# I wonder if we can define a single maneuver that takes us on a path from Earth all the way to Saturn.
# I imagine we could specify target orbits rather than simply target points based on ephemeris files,
# but that doesn't work for the next mission.
jupiter_end, jupiter_cost, jupiter_ephem = assist_to_planet(
    Jupiter, jupiter_date, init_orbit_end, "black", (0, 150)
)
saturn_end, saturn_cost, saturn_ephem = assist_to_planet(
    Saturn, saturn_date, jupiter_end, "blue", (30, 80)
)
uranus_end, uranus_cost, uranus_ephem = assist_to_planet(
    Uranus, uranus_date, saturn_end, "red", (5, 70)
)
neptune_end, neptune_cost, neptune_ephem = assist_to_planet(
    Neptune, neptune_date, uranus_end, "green", (40, 68)
)

vectors = [
    magnitude(init_orbit),
    magnitude(jupiter_end),
    magnitude(saturn_end),
    magnitude(uranus_end),
    magnitude(neptune_end),
]

pprint(vectors)
pprint(f"Jupiter: {jupiter_cost}")
pprint(f"Saturn: {saturn_cost}")
pprint(f"Uranus: {uranus_cost}")
pprint(f"Neptune: {neptune_cost}")
