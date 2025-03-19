#! venv/bin/python3

import numpy as np
import scipy as sp

import yaml
import warnings

from mpi4py import MPI
from mpi4py.util import dtlib


import matplotlib.pyplot as plt
from myplot import *


eps = 1e-8
max_iter = 100


class Primary:

    def __init__(self, a, e):
        self.e = e
        self.a = a

        if e < 0 or e > 1:
            raise AttributeError("Eccentricity must lie within 0 to 1 range for a closed orbit")


    def solve_kepler(self, M):
        e = self.e
        E0 = 0
        E = M

        i = 0

        while np.mean(np.abs(E - E0)) > 1e-8:
            E0 = E
            E = e * np.sin(E0) + M 

            i += 1

            if i > max_iter:
                warnings.warn(f"Warning! Kepler solver has not converged after {max_iter} iterations")

        return E % (2*np.pi)


    def __call__(self, t):
        """
        A function to calculate the distance between system's 
        center of mass and primary body at moments t 
        (given semi-major orbit axis a and eccentrisity e)

        :param t: float scalar of array-like moments of time at which evaluate primary distance to the system's center of mass
        :return r: array-like evaluated distances 
        """
        M = t   # according to M = 2\pi t/P, and P = 2\pi
        E = self.solve_kepler(M)

        r = self.a * (1.0 - self.e * np.cos(E)) / 2

        return r

class DE:

    def __init__(self, r, track_progress = False):
        self.r = r
        self.track_progress = track_progress

    def ivp(self, t, y):
        deriv = np.empty_like(y)
        
        z = y[0::2]
        v = y[1::2]

        deriv[0::2] = v
        deriv[1::2] = -z / np.sqrt(self.r(t)**2 + z**2) ** 3

        if self.track_progress and size > 1:
            comm.isend(t, dest = 0, tag=tag_data)

        return deriv
        

def update_progress(thread, t):
    nlines = size - thread

    print(f"\033[{nlines}A", end="")
    print(f"\rTHREAD {thread}: {t/tf*100:.1f} %")

    if nlines > 1:
        print(f"\033[{nlines-1}B", end="")



def run_parallel(y0, rank, size):
    # find how much realizations are
    # expected to be calculated by a thread
    num_per_thread = int(nz * ndot / size)


    # split initial conditions to threads
    # add additional realizations to the last thread
    # if there are leftover ones
    if rank < size:
        y0 = y0[2*(rank-1) * num_per_thread :\
                2*rank * num_per_thread]
    else:
        y0 = y0[2*(rank-1) * num_per_thread:]

    
    # perform the integration
    res = sp.integrate.solve_ivp(de.ivp, (t0, tf), y0, 
                                 method = "Radau",
                                 t_eval = t_eval,
                                 )

                                 #max_step = 1e-1)

    comm.send(None, dest = 0, tag=tag_end)


    # split the result for convinience
    t = res.t
    zs = res.y[0::2, :]
    vs = res.y[1::2, :]

    # save the result into according files
    np.savetxt(f"{rank}t.dat", t)
    np.savetxt(f"{rank}zs.dat", zs)
    np.savetxt(f"{rank}vs.dat", vs)


if __name__ == "__main__":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    datatype = MPI.FLOAT
    itemsize = datatype.Get_size()

    win_size = itemsize if rank == 0 else 0
    win = MPI.Win.Allocate(win_size, comm=comm)
    np_dtype = dtlib.to_numpy_dtype(datatype)
    buf = np.empty(1, dtype=np_dtype)

    tag_data = 42
    tag_end = 23



    with open("config.yaml") as stream:
        conf = yaml.safe_load(stream)

    conf["parallel"].update({"size" : size})

    if size == 1:
        track = False

    # save used configuration to use it 
    # while plotting the results
    if rank == 0 or size == 1:
        with open(".last_config.yaml", "w") as stream:
            yaml.safe_dump(conf, stream)

    # read some configuration parameters
    a = conf["orbit"]["a"]
    e = conf["orbit"]["e"]
    t0 = conf["grid"]["t0"] * 2*np.pi
    tf = conf["grid"]["tf"] * 2*np.pi 
    z0 = conf["grid"]["z0"]
    z1 = conf["grid"]["z1"]
    nz = conf["grid"]["nz"]
    dotz0 = conf["grid"]["dotz0"]
    dotz1 = conf["grid"]["dotz1"]
    ndot = conf["grid"]["ndot"]
    track = conf["parallel"]["track"]

    if t0 < 0 or tf < 0:
        raise AttributeError("Time domain can only be positive")



    # make a grid of initial conditions
    z = np.linspace(z0, z1, nz, dtype = np.float128)
    dotz = np.linspace(dotz0, dotz1, ndot, dtype = np.float128)

    t_eval = np.arange(t0, tf, 2*np.pi)

    z, dotz = np.meshgrid(z, dotz)

    # merge the grid into 1d-array
    # structured as:
    # (z0, v0, z1, v1, . . . , zN, vN)
    y0 = np.empty(2*nz*ndot)
    y0[0::2] = z.flatten()
    y0[1::2] = dotz.flatten()


    # declare primary body and the differential equation
    body = Primary(a, e)
    de = DE(body, track_progress = track)

    # solve in parallel
    if size > 1:

        # update the progress using 0th thread
        if rank == 0 and track:
            remaining = comm.Get_size() - 1

            print(f"THREAD {0}: tracking progress")
            for i in range(size - 1):
                print(f"THREAD {i+1}: 0.0 %")

            while remaining > 0:
                s = MPI.Status()
                comm.Probe(status = s)

                if s.tag == tag_data:
                    t = comm.recv(source = s.source, tag = tag_data)
                    update_progress(s.source, t)

                elif s.tag == tag_end:
                    comm.recv(source = s.source, tag = tag_end)
                    remaining -= 1

        # and solve the problem with other threads
        elif rank > 0 and track:
            run_parallel(y0, rank, size-1)

        # use all threads for calculation
        elif not track:
            run_parallel(y0, rank+1, size)

    # solve on a single core
    else:
        run_parallel(y0, 1, 1)



