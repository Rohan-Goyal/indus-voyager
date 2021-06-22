# Two lines describe a particular time
import numpy as np
from datetime import datetime as dt


class datapoint:
    def __init__(self, epoch, v):
        self.epoch = epoch
        self.v = v

    def __repr__(self):
        return str([self.epoch.isoformat(), self.v])

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
            text = text.replace(f"V{i}=", "")
        return np.array([float(i) * 1000 for i in text.split()])

    @classmethod
    def from_lines(self, line1, line2):
        epoch = self.date_cleaned(line1)
        v_list = self.v_cleaned(line2)
        return datapoint(epoch, v_list)


def _from_file(fname):
    # Returns a list of datapoints
    with open(fname) as filename:
        lines = filename.readlines()
    start = lines.index("$$SOE\n") + 1
    end = lines.index("$$EOE\n") - 1
    relevant = lines[start:end]
    for i in range(0, len(relevant) - 1, 2):
        yield datapoint.from_lines(relevant[i].strip(), relevant[i + 1].strip())


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

    # I want something that basically tells us if change is > or < than average.


# ds=from_file("/Users/rohangoyal/Downloads/may25.txt")
ds = from_file("/Users/rohangoyal/Downloads/horizons_results (10).txt")
