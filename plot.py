#! venv/bin/python3

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle

from myplot import *

import yaml

set_tex()

if __name__ == "__main__":

    with open(".last_config.yaml", "r") as stream:
        conf = yaml.safe_load(stream)

    size = conf["parallel"]["size"]
    rev_min = conf["plotter"]["rev_min"]
    track = conf["parallel"]["track"]



    plt.figure(figsize = (6, 6))
    plt.title("Phase portrait")
    plt.xlabel(r"$z$")
    plt.ylabel(r"$\dot z$")

    plt.xlim(-2.5, 2.5)
    plt.ylim(-2, 2)

    if size > 1:
        size -= 1

    for i in range(size):
        print(f"\rProcessing thread {i+1}", end = "")

        ts = np.loadtxt(f"{i+1}t.dat")
        zs = np.loadtxt(f"{i+1}zs.dat")
        vs = np.loadtxt(f"{i+1}vs.dat")

        zs = np.atleast_2d(zs)
        vs = np.atleast_2d(vs)

        n_bodies, n_times = zs.shape
        bodies = np.arange(n_bodies)

        i = np.where(np.logical_and(
                np.mod(ts, 2*np.pi) < 1e-2,
                ts > rev_min * 2*np.pi))[0]

#        ts_new = ts[i]
#        zs_new = zs[:, i]
#        vs_new = vs[:, i]


#        ts_new = np.arange(rev_min, ts.max()-1, 1)
#
#        zs_new = sp.interpolate.interpn(points = (bodies, ts),
#                                      values = zs, 
#                                      xi = np.meshgrid(bodies, ts_new),
#                                      method = "cubic").T
#
#        vs_new = sp.interpolate.interpn(points = (bodies, ts),
#                                      values = vs, 
#                                      xi = np.meshgrid(bodies, ts_new),
#                                      method = "cubic").T
#

        #plt.plot(zs_new, vs_new,
        plt.plot(zs, vs,
                 ",",
                 color = "black", alpha = 1,
                 )
        
    print("")

    save_image("123.pdf")





