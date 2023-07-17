"""
Simple demonstration of the Pool() class within multiprocessing.
"""
import os
import time
import shutil
import numpy as np
import subprocess as sp
import multiprocessing as mp
import matplotlib.pyplot as plt


def fac(n: int) -> int:
    """Return factorial."""
    fac = 1
    for i in range(n, 1, -1):
        fac *= i
    return fac


def Mq(q: int) -> int:
    """Factorial term in Chudnovsky."""
    return fac(6*q)//(fac(3*q)*fac(q)**3)


def Lq(q: int) -> int:
    """Linear term in Chudnovsky."""
    return 545140134 * q + 13591409


def Xq(q: int) -> int:
    """Exponential term in Chudnovsky."""
    return int(-262537412640768000)**q


def piChudnovsky(n: int) -> None:
    """Chudnovsky algorithm for computing pi."""
    print(f"Chudnovsky with n={n}")
    C = 426880 * np.sqrt(10005)
    pi = 0
    for i in range(n):
        pi += (Mq(i)*Lq(i)/Xq(i))
    print(f"Chudnovsky with n={n}... Done! pi={C/pi}")


if __name__ == "__main__":

    # Monitor with btop
    print("Monitoring with btop in separate window!\n\n")
    sp.run(["xfce4-terminal " +
            "--geometry=150x60-0+0 " +
            "--zoom=-2 " +
            "-e 'btop' " +
            ">/dev/null 2>/dev/null &"],
           shell=True, check=False, capture_output=True)
    time.sleep(5)

    # Example with computation in Python
    print("Python example: Chudnovsky algorithm for computing pi")
    with mp.Pool(processes=8) as p:
        p.map(piChudnovsky, [i*100 for i in range(10, 25)])
        p.close()
        p.join()
    print("\n\n")

    # Example with an external executable
    def runFortran(run: int) -> None:
        """Set up and run embarrassingly parallel FORTRAN code."""
        os.mkdir(f"runs/run{run:0>2}")
        shutil.copy("fmap/fmap.x", f"runs/run{run:0>2}/")
        shutil.copy("fmap/input", f"runs/run{run:0>2}/")
        print(f"Running fmap.x in runs/run{run:0>2}... ")

        sp.run([f"./fmap.x > run{run:0>2}.log"],
               shell=True, cwd=f"runs/run{run:0>2}")

        print(f"Running fmap.x in runs/run{run:0>2}... Done")

    print("Executable example: Semiclassical dynamics with FORTRAN")
    os.mkdir("runs")
    with mp.Pool(processes=8) as p:
        p.map(runFortran, [i+1 for i in range(16)])
        p.close()
        p.join()
    dat = np.genfromtxt("runs/run01/Cpop.out")
    for i in range(1, 16):
        dat += np.genfromtxt(f"runs/run{i+1:0>2}/Cpop.out")
    dat = dat / 16

    plt.plot(dat[1:, 0], dat[1:, 1])
    plt.plot(dat[1:, 0], dat[1:, 2])
    plt.show()
