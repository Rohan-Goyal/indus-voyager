import pygmo as pg
from pykep import epoch, epoch_from_string
from pykep.planet import jpl_lp
from pykep.trajopt import mga_1dsm, mga


y=365.25

def ex5():

    # We define an Earth-Venus-Earth problem (single-objective)
    seq = [jpl_lp('earth'), jpl_lp('jupiter'), jpl_lp('saturn'), jpl_lp('uranus'),jpl_lp('neptune')]
    udp = mga_1dsm(
        seq=seq,
        t0=[epoch_from_string("1977-08-20 14:25:00"), epoch_from_string("1977-08-20 14:35:00")],
        tof=[[1.5*y,2 * y], [1.8*y,2.2* y], [4.3*y,4.8*y], [3.5*y,3.8*y]], # TODO: Lower bounds
        vinf=[35,40],
        add_vinf_dep=False,
        add_vinf_arr=False,
        multi_objective=False,
        tof_encoding="direct"
    )

    pg.problem(udp)
    # We solve it!!
    uda = pg.sade(gen=100)
    archi = pg.archipelago(algo=uda, prob=udp, n=8, pop_size=20)
    print(
        "Running a Self-Adaptive Differential Evolution Algorithm .... on 8 parallel islands")
    archi.evolve(10)
    archi.wait()
    sols = archi.get_champions_f()
    idx = sols.index(min(sols))
    print("Done!! Solutions found are: ", archi.get_champions_f())
    udp.pretty(archi.get_champions_x()[idx])
    udp.plot(archi.get_champions_x()[idx])

def run_sim():
    seq = [jpl_lp('earth'), jpl_lp('jupiter'), jpl_lp('saturn'), jpl_lp('uranus'),jpl_lp('neptune')]
    udp = mga(
        seq=seq,
        t0=[epoch_from_string("1977-08-20 14:25:00"), epoch_from_string("1977-08-20 14:35:00")],
        tof=[[1.5*y,2 * y], [1.8*y,2.2* y], [4.3*y,4.8*y], [3.5*y,3.8*y]], # TODO: Lower bounds
        vinf=40.0,
        tof_encoding="direct"
    )

    pg.problem(udp)
    # We solve it!!
    uda = pg.sade(gen=100)
    archi = pg.archipelago(algo=uda, prob=udp, n=10, pop_size=25)
    print(
        "Running a Self-Adaptive Differential Evolution Algorithm .... on 8 parallel islands")
    archi.evolve(15)
    archi.wait()
    sols = archi.get_champions_f()
    idx = sols.index(min(sols))
    print("Done!! Solutions found are: ", archi.get_champions_f())
    udp.pretty(archi.get_champions_x()[idx])
    udp.plot(archi.get_champions_x()[idx])


def modern():

    # We define an Earth-Venus-Earth problem (single-objective)
    seq = [jpl_lp('earth'), jpl_lp('jupiter'), jpl_lp('uranus'),jpl_lp('neptune')]
    udp = mga_1dsm(
        seq=seq,
        t0=[epoch_from_string("2024-01-01 00:00:00"), epoch_from_string("2025-12-31 23:59:59")],
        tof=[[y,3*y], [4*y,8*y], [2.5*y,4.5*y]], # TODO: Lower bounds
        vinf=[35,40],
        add_vinf_dep=False,
        add_vinf_arr=False,
        multi_objective=False,
        tof_encoding="direct"
    )

    pg.problem(udp)
    # We solve it!!
    uda = pg.sade(gen=120)
    archi = pg.archipelago(algo=uda, prob=udp, n=10, pop_size=25)
    print(
        "Running a Self-Adaptive Differential Evolution Algorithm .... on 8 parallel islands")
    archi.evolve(12)
    archi.wait()
    sols = archi.get_champions_f()
    idx = sols.index(min(sols))
    print("Done!! Solutions found are: ", archi.get_champions_f())
    udp.pretty(archi.get_champions_x()[idx])
    udp.plot(archi.get_champions_x()[idx])

def modern_unpowered():

    # We define an Earth-Venus-Earth problem (single-objective)
    seq = [jpl_lp('earth'), jpl_lp('jupiter'), jpl_lp('uranus'),jpl_lp('neptune')]
    udp = mga(
        seq=seq,
        t0=[epoch_from_string("2024-01-01 00:00:00"), epoch_from_string("2025-12-31 23:59:59")],
        tof=[[y,3*y], [4*y,8*y], [2.5*y,4.5*y]], # TODO: Lower bounds
        vinf=[35,40],
    )

    pg.problem(udp)
    # We solve it!!
    uda = pg.sade(gen=120)
    archi = pg.archipelago(algo=uda, prob=udp, n=8, pop_size=20)
    print(
        "Running a Self-Adaptive Differential Evolution Algorithm .... on 8 parallel islands")
    archi.evolve(12)
    archi.wait()
    sols = archi.get_champions_f()
    idx = sols.index(min(sols))
    print("Done!! Solutions found are: ", archi.get_champions_f())
    udp.pretty(archi.get_champions_x()[idx])
    udp.plot(archi.get_champions_x()[idx])
if __name__=="__main__":
    modern_unpowered()
