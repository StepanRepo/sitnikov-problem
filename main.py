#! venv/bin/python3

import numpy as np
import scipy as sp

import yaml
import warnings

from mpi4py import MPI
from mpi4py.util import dtlib


import matplotlib.pyplot as plt


eps = 1e-8
max_iter = 100


class Primary:
    """
    A class representing a primary body in the Sitnikov problem.

    Parameters
    ----------
    a : float
        Semi-major axis of the orbit
    e : float
        Eccentricity of the orbit (must be between 0 and 1)
    """

    def __init__(self, a, e):
        self.e = e
        self.a = a

        if e < 0 or e > 1:
            raise AttributeError("Eccentricity must lie within 0 to 1 range for a closed orbit")


    def solve_kepler(self, M):
        """
        Solves Kepler's equation iteratively for eccentric anomaly.

        Parameters
        ----------
        M : float or numpy.ndarray
            Mean anomaly

        Returns
        -------
        float or numpy.ndarray
            Eccentric anomaly solution modulo 2π
        """
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
        Calculate the distance between system's center of mass and primary body.

        Parameters:
        -----------
        t : float or array-like
            Moments of time at which to evaluate primary distance

        Returns:
        --------
        float or array-like
            Distances between the system's center of mass and the primary body
        """
        M = t   # according to M = 2\pi t/P, and P = 2\pi
        E = self.solve_kepler(M)

        r = self.a * (1.0 - self.e * np.cos(E)) / 2

        return r

class DE:
    """
    Differential equation solver for the Sitnikov problem.

    Implements the equations of motion for a third body in the Sitnikov problem,
    where two primary bodies move in elliptical orbits.
    """

    def __init__(self, r, track_progress = False):
        """
        Initialize the differential equation solver.

        Parameters:
        -----------
        r : callable
            Function that returns the distance of primary bodies from center of mass
        track_progress : bool, optional
            Whether to track and report integration progress
        """
        self.r = r
        self.track_progress = track_progress

    def ivp(self, t, y):
        """
        Define the initial value problem for the Sitnikov system.

        Implements the equations of motion for the third body in the Sitnikov problem.

        Parameters:
        -----------
        t : float
            Current time value
        y : array-like
            Current state vector [z1, v1, z2, v2, ..., zn, vn]

        Returns:
        --------
        array-like
            Derivatives of state vector [v1, a1, v2, a2, ..., vn, an]
        """
        deriv = np.empty_like(y)
        
        z = y[0::2]
        v = y[1::2]

        deriv[0::2] = v
        deriv[1::2] = -z / np.sqrt(self.r(t)**2 + z**2) ** 3

        if self.track_progress and size > 1:
            if np.mod(t, 1) < 1e-2:
                comm.isend(t, dest = 0, tag=tag_data)

        return deriv
        

def update_progress(thread, t):
    """
    Update and display integration progress for a specific thread.

    Parameters:
    -----------
    thread : int
        Thread ID for which progress is being updated
    t : float
        Current time in the integration
    """

    nlines = size - thread

    print(f"\033[{nlines}A", end="")
    print(f"\rTHREAD {thread}: {t/tf*100:.1f} %")

    if nlines > 1:
        print(f"\033[{nlines-1}B", end="")



def run_parallel(y0, rank, size, *args, **kwargs):
    """
    Run the integration in parallel using multiple threads.

    Splits initial conditions among threads, performs integration,
    and saves results to files.

    Parameters:
    -----------
    y0 : array-like
        Initial conditions for all integrations
    rank : int
        Rank of the current MPI process
    size : int
        Total number of processes for computation
    *args, **kwargs
        Additional arguments passed to scipy.integrate.solve_ivp
    """

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
                                 *args, **kwargs)

    if track:
        comm.send(None, dest = 0, tag=tag_end)


    # split the result for convinience
    t = res.t
    zs = res.y[0::2, :]
    vs = res.y[1::2, :]

    # save the result into according files
    np.savetxt(f"result/{rank}t.dat", t)
    np.savetxt(f"result/{rank}zs.dat", zs)
    np.savetxt(f"result/{rank}vs.dat", vs)


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
    method = conf["integrator"]["method"]

    if t0 < 0 or tf < 0:
        raise AttributeError("Time domain can only be positive")



    # make a grid of initial conditions
    z = np.linspace(z0, z1, nz, dtype = np.float64)
    dotz = np.linspace(dotz0, dotz1, ndot, dtype = np.float64)

    # make a time grid to evaluate function at 
    t_eval = np.arange(t0, tf, 2*np.pi)

    # make a grid of initial parameters
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

    # solve on a single core
    if size == 1:
        print(f"Running in a single thread without tracking")
        run_parallel(y0, 1, 1, method = method, t_eval = t_eval)
    # solve in parallel
    else:

        # update the progress using 0th thread
        if rank == 0 and track:
            print(f"Running in {size} threads with tracking")

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
            run_parallel(y0, rank, size-1, 
                         method = method,
                         t_eval = t_eval)

        # use all threads for calculation
        elif not track:
            print(f"Running in {size} threads without tracking")
            run_parallel(y0, rank+1, size,
                         method = method,
                         t_eval = t_eval)




