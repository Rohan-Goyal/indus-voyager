# Two lines describe a particular time
import numpy as np
from datetime import datetime as dt
from astropy import units as u

class datapoint:
    def __init__(self, epoch, r, v):
        self.epoch = epoch
        self.r=r
        self.v = v

    def __repr__(self):
        return str([self.epoch.isoformat(), self.r, self.v])

    def date_cleaned(text):
        val = text.split("A.D. ")[1].replace(" TDB", "").replace("0000", "000")
        months = [
            "Jan",
            "Feb",
            "Mar",
            "Apr",
            "May",
            "Jun",
            "Jul",
            "Aug",
            "Sep",
            "Oct",
            "Nov",
            "Dec",
        ]
        for i, v in enumerate(months):
            val = (
                val.replace(v, "0" + str(i + 1))
                if i < 9
                else val.replace(v, str(i + 1))
            )
        return dt.fromisoformat(val)

    def v_cleaned(text):
        for i in "XYZ":
            text = text.replace(f"V{i}=", "").replace(f"{i} =","")
        return np.array([float(i) * 1000 for i in text.split()])

    @classmethod
    def from_lines(self, line1, line2, line3):
        epoch = self.date_cleaned(line1)
        r_list=self.v_cleaned(line2)
        v_list = self.v_cleaned(line3)
        return datapoint(epoch, r_list, v_list)


def _from_file(fname):
    # Returns a list of datapoints
    with open(fname) as filename:
        lines = filename.readlines()
    start = lines.index("$$SOE\n") + 1
    end = lines.index("$$EOE\n") - 2
    relevant = lines[start:end]
    for i in range(0, len(relevant) - 1, 3):
        yield datapoint.from_lines(relevant[i].strip(), relevant[i + 1].strip(), relevant[i+2].strip())


def from_file(fname):
    return list(_from_file(fname))


# Comparing datapoints:
# Method to take magnitude of velocity
# Method to do vector subtraction
def _vsub(v1, v2):
    for i in range(len(v1)):
        yield v2[i] - v1[i]


def vsub(v1, v2):
    return list(_vsub(v1, v2))


def norm(vec):
    return np.sqrt(sum([i ** 2 for i in vec]))


def vdiff(v1, v2):
    return norm(vsub(v1, v2))


def maxdiff(dataset):
    m = 0
    m_idx = -1
    for i in range(len(dataset) - 1):
        diff = vdiff(dataset[i].v, dataset[i + 1].v)
        if abs(diff) >= abs(m):
            m = diff
            m_idx = i
    return m, m_idx


def mean_dv(data):
    return vdiff(data[0].v, data[-1].v) / len(data)


def above(data):
    return [
        data[i]
        for i in range(1, len(data))
        if vdiff(data[i - 1].v, data[i].v) > 1.02 * mean_dv(data)
    ]


def vdiffs(data):
    diffdata = [0] * (len(data) - 1)
    for i in range(len(data) - 1):
        diff = vdiff(data[i].v, data[i + 1].v)
        diffdata[i] = (data[i].epoch, diff)
    return diffdata


def vnorms(data):
    normdata = [0] * len(data)
    for i in range(len(data)):
        diff = norm(data[i].v)
        normdata[i] = (data[i].epoch, diff)
    return normdata

def rnorms(data):
    normdata = [0] * len(data)
    for i in range(len(data)):
        diff = norm(data[i].r)
        normdata[i] = (data[i].epoch, diff)
    return normdata
    # I want something that basically tells us if change is > or < than average.

def closest_r(data, const): # Find epoch/rv of object which is closest to the position vector const
    diffs=[vdiff(const, i.r) for i in data]
    idx=diffs.index(min(diffs))
    return data[idx]

def closest_v(data,const):
    diffs=[vdiff(const, i.v) for i in data]
    idx=diffs.index(min(diffs))
    return data[idx]

def closest_rv(data, constr, constv):
    rdiffs=[vdiff(constr, i.r) for i in data]
    vdiffs=[vdiff(constv, i.v) for i in data]
    diffs=[rdiffs[i]+vdiffs[i] for i in range(len(data))]
    idx=diffs.index(min(diffs))
    return data[idx]

# ds=from_file("/Users/rohangoyal/Downloads/may25.txt")
# Jupiter July 1 data:
# 2443325.500000000 = A.D. 1977-Jul-01 00:00:00.0000 TDB
#  X = 1.551783168059400E+08 Y = 8.877133948602245E+08 Z =-7.056679235727489E+06
#  VX=-4.175487351667007E+01 VY=-1.028913798527010E+00 VZ= 2.726943920685317E-01

voy2_jup=datapoint.from_lines("2443325.500000000 = A.D. 1977-Jul-01 00:00:00.0000 TDB",
                              "X = 1.551783168059400E+08 Y = 8.877133948602245E+08 Z =-7.056679235727489E+06",
                              "VX=-4.175487351667007E+01 VY=-1.028913798527010E+00 VZ= 2.726943920685317E-01"
                              )

# Uranus Jan 21 data:
# 2446451.500000000 = A.D. 1986-Jan-21 00:00:00.0000 TDB
# X =-4.754062367113435E+08 Y =-2.933414590935786E+09 Z =-3.268767410831213E+06
# VX= 3.295603995565554E+01 VY= 1.347721156523741E+01 VZ= 3.656038905200099E-02

voy2_ur= datapoint.from_lines("2446451.500000000 = A.D. 1986-Jan-21 00:00:00.0000 TDB",
                              "X =-4.754062367113435E+08 Y =-2.933414590935786E+09 Z =-3.268767410831213E+06",
                              "VX= 3.295603995565554E+01 VY= 1.347721156523741E+01 VZ= 3.656038905200099E-02"
                              )

ds_jup=from_file("/Users/rohangoyal/Downloads/jupiter 2024.txt")
ds_ur=from_file("/Users/rohangoyal/Downloads/uranus 2034.txt")

# In [1]: closest_rv(ds,voy2_jup.r, voy2_jup.v)
# Out[1]: ['2026-04-01T00:00:00', array([-2.01609233e+11,  7.31228899e+11,  4.87158347e+09]), array([-17026.5537969 ,  24154.0927399 ,    284.36336515])]
