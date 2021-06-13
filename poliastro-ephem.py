from pprint import pprint
from astropy import units as u
from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time, TimeDelta
from poliastro.bodies import Earth, Jupiter, Neptune, Saturn, Sun, Uranus
from poliastro.ephem import Ephem
from poliastro.plotting import StaticOrbitPlotter
from poliastro.twobody import Orbit
from poliastro.util import time_range, norm

solar_system_ephemeris.set("jpl")

full_time = time_range("1977-08-21", end="1989-08-30", periods=300)
voy_ephem = Ephem.from_horizons("Voyager 2", full_time, id_type="majorbody")

vectors = voy_ephem.rv(epochs=full_time)

table = [[i for i in full_time]]
table.append(vectors[0])
table.append([i.to(u.km / u.s) for i in vectors[1]])
table.append([norm(i) for i in vectors[0]])
table.append([norm(i).to(u.km / u.s) for i in vectors[1]])

# Each column (table[0,i]], table[1,i], table[2,i]) gives us an epoch,position,velocity tuple.


def column(i):
    for x in range(0, len(table)):
        yield table[x][i]


def col(i):
    # pprint(list(column(i)))
    return list(column(i))


def vsub(i, j):
    if i >= j:
        return col(i)[2] - col(j)[2], norm(col(i)[2] - col(j)[2])
    return col(j)[2] - col(i)[2], norm(col(j)[2] - col(i)[2])