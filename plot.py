#! venv/bin/python3

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle

from myplot import *

import yaml
from glob import glob


if __name__ == "__main__":
    try:
        plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.sans-serif": "serif",
        "font.size"   : 12
        })
    except:
        print("TeX font is not availible")




    with open("config.yaml", "r") as stream:
        conf = yaml.safe_load(stream)

    rev_min = conf["plotter"]["rev_min"]
    xmin  = conf["plotter"]["zmin"]
    xmax  = conf["plotter"]["zmax"]
    ymin  = conf["plotter"]["dotzmin"]
    ymax  = conf["plotter"]["dotzmax"]
    dpi  = conf["plotter"]["dpi"]
    track = conf["parallel"]["track"]

    res_list = glob("result/*t.dat")
    size = len(res_list)





    plt.figure(figsize = (6, 6))
    plt.title("Phase portrait")
    plt.xlabel(r"$z$")
    plt.ylabel(r"$\dot z$")

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)


    for i in range(size):
        print(f"\rProcessing thread {i+1}", end = "")

        ts = np.loadtxt(f"result/{i+1}t.dat")
        zs = np.loadtxt(f"result/{i+1}zs.dat")
        vs = np.loadtxt(f"result/{i+1}vs.dat")

        zs = np.atleast_2d(zs)
        vs = np.atleast_2d(vs)

        plt.plot(zs, vs,
                 ".",
                 color = "black"
                 )
        
    print("")

    plt.savefig("123.png", dpi = dpi)






