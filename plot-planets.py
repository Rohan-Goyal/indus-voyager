import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import solar_system_ephemeris

solar_system_ephemeris.set("jpl")
from poliastro.bodies import Sun, Earth, Mars, Jupiter, Saturn, Uranus, Neptune
from poliastro.plotting import StaticOrbitPlotter
from poliastro.util import norm, time_range

launch = Time("2024-01-01 00:00", scale="utc").tdb
plotter = StaticOrbitPlotter()

# def orbit(planet, duration, start=launch):
#     return Ephem.from_body(planet, time_range(start, start+duration))


def years(n):
    return TimeDelta(86400.0 * 365 * n, format="sec")


def months(n):
    return TimeDelta(86400.0 * 30 * n, format="sec")


def plot_planet(planet, duration, start=launch, plotter=plotter):
    plotter.plot_body_orbit(
        planet, start + years(duration), label=f"{planet} at time of flyby"
    )


def plot_all(start=launch, plotter=plotter):
    plot_planet(Earth, 0, start, plotter)
    plot_planet(Mars, 1, start, plotter)
    plot_planet(Jupiter, 2, start, plotter)
    plot_planet(Saturn, 4, start, plotter)
    plot_planet(Uranus, 9.5, start, plotter)
    plot_planet(Neptune, 13, start, plotter)


def multiple(start=launch, inc=3, it=8):
    for i in range(0, it):
        start += months(inc)
        plot_all(start, StaticOrbitPlotter())


multiple()


def plot_naive(start=launch, plotter=StaticOrbitPlotter()):
    plot_planet(Earth, 0, start)
    plot_planet(Jupiter, 0, start)
    plot_planet(Saturn, 0, start)
    plot_planet(Uranus, 0, start)
    plot_planet(Neptune, 0, start)
