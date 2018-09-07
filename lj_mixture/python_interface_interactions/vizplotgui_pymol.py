#!/usr/bin/env python -i
# preceeding line should have path for Python on your machine

# vizplotgui_pymol.py
# Purpose: viz running LAMMPS simulation via PyMol with plot and GUI
# Syntax:  vizplotgui_pymol.py in.lammps Nfreq compute-ID
#          in.lammps = LAMMPS input script
#          Nfreq = plot data point and viz shapshot every this many steps
#          compute-ID = ID of compute that calculates temperature
#                       (or any other scalar quantity)


from __future__ import print_function

import os
import sys
import time

try:
    path = os.environ["LAMMPS_PYTHON_TOOLS"]
except KeyError:
    path = "../../libs/pizza"
sys.path.append(path)

try:
    PYMOL_PATH = os.environ['PYMOL_PATH']
except KeyError:
    PYMOL_PATH = '/usr/lib/python2.7/dist-packages/pymol'
    os.environ['PYMOL_PATH'] = PYMOL_PATH


# methods called by GUI

def run():
    global runflag
    runflag = 1


def stop():
    global runflag
    runflag = 0


def set_temp(value):
    global target_temp
    target_temp = slider.get()


def quit():
    global breakflag
    breakflag = 1

def set_eps11(value):
  global target_eps11
  target_eps11 = slider_eps11.get()


def set_eps12(value):
  global target_eps12
  target_eps12 = slider_eps12.get()


def set_eps22(value):
  global target_eps22
  target_eps22 = slider_eps22.get()

# method called by timestep loop every Nfreq steps
# read dump snapshot and viz it, update plot with compute value

def update(ntimestep):
    d.next()
    d.unscale()
    p.single(ntimestep)
    pm.load("tmp.pdb")
    pm.forward()
    value = lmp.extract_compute(compute, 0, 0)
    xaxis.append(ntimestep)
    yaxis.append(value)
    gn.plot(xaxis, yaxis)


# parse command line

argv = sys.argv
if len(argv) != 3:
    print("Syntax: vizplotgui_pymol.py in.lammps Nfreq")
    sys.exit()

infile = sys.argv[1]
nfreq = int(sys.argv[2])
compute = "thermo_temp"

me = 0

from mpi4py import MPI

comm = MPI.COMM_WORLD
me = comm.Get_rank()
nprocs = comm.Get_size()

from lammps import lammps

lmp = lammps()

# run infile all at once
# assumed to have no run command in it
# dump a file in native LAMMPS dump format for Pizza.py dump tool

lmp.file(infile)
lmp.command("thermo %d" % nfreq)
lmp.command("dump python all atom %d tmp.dump" % nfreq)

# initial 0-step run to generate initial 1-point plot, dump file, and image

lmp.command("run 0 pre yes post no")
value = lmp.extract_compute(compute, 0, 0)
ntimestep = 0
xaxis = [ntimestep]
yaxis = [value]

breakflag = 0
runflag = 0
target_temp = value

# define interactions (epsilon)
target_eps11 = 1.0
target_eps12 = 1.0
target_eps22 = 1.0

str_pair_coeff_11 = "pair_coeff	1 1 %g 1 2.5"
str_pair_coeff_12 = "pair_coeff	1 2 %g 1 2.5"
str_pair_coeff_22 = "pair_coeff 2 2 %g 1 2.5"

# wrapper on PyMol
# just proc 0 handles reading of dump file and viz

if me == 0:
    import pymol

    pymol.finish_launching()

    from dump import dump
    from pdbfile import pdbfile
    from pymol import cmd as pm

    d = dump("tmp.dump", 0)
    p = pdbfile(d)
    d.next()
    d.unscale()
    p.single(ntimestep)
    pm.load("tmp.pdb")
    pm.show("spheres", "tmp")
    pm.color("red", "name 1")
    pm.select("type_1", "name 1", enable=0)
    pm.color("blue", "name 2")
    pm.select("type_2", "name 2", enable=0)

# display GUI with run/stop buttons and slider for temperature

if me == 0:
    # fix Tkinter import error
    try:
        from Tkinter import *
    except:
        from tkinter import *

    tkroot = Tk()
    tkroot.withdraw()
    root = Toplevel(tkroot)
    root.title("LAMMPS GUI")

    frame = Frame(root)
    Button(frame, text="Run", command=run).pack(side=LEFT)
    Button(frame, text="Stop", command=stop).pack(side=LEFT)

    # slider for temperature
    slider = Scale(frame, from_=0.0, to=3.0, resolution=0.1,
                   orient=HORIZONTAL, label="Temperature")
    slider.bind('<ButtonRelease-1>', set_temp)
    slider.set(target_temp)
    slider.pack(side=LEFT)
    Button(frame, text="Quit", command=quit).pack(side=RIGHT)
    frame.pack()

    label = Label(root, text="Interactions")
    label.pack(fill=X, pady=10, padx=10)

    frame = Frame(root)
    # slider for epsilon 1 - 1
    Label(frame, text="Epsilon red - red").grid(row=0,column=0)
    value_eps11 = IntVar()
    slider_eps11 = Scale(frame, from_=0.0, to=2.0, resolution=0.05,
                         orient=HORIZONTAL, variable=value_eps11, showvalue=0)
    slider_eps11.bind('<ButtonRelease-1>', set_eps11)
    slider_eps11.set(target_eps11)
    slider_eps11.grid(row=0,column=1)
    Label(frame, textvariable=value_eps11).grid(row=0, column=2)

    # # slider for epsilon 2 - 2
    Label(frame, text="Epsilon blue - blue").grid(row=1, column=0)
    value_eps22 = IntVar()
    slider_eps22 = Scale(frame, from_=0.0, to=2.0, resolution=0.05,
                          orient=HORIZONTAL, variable=value_eps22, showvalue=0)
    slider_eps22.bind('<ButtonRelease-1>', set_eps22)
    slider_eps22.set(target_eps22)
    slider_eps22.grid(row=1,column=1)
    Label(frame, textvariable=value_eps22).grid(row=1, column=2)

    # slider for epsilon 1 - 2
    Label(frame, text="Epsilon red - blue").grid(row=2, column=0)
    value_eps12 = IntVar()
    slider_eps12 = Scale(frame, from_=0.0, to=2.0, resolution=0.05,
                          orient=HORIZONTAL, variable=value_eps12, showvalue=0)
    slider_eps12.bind('<ButtonRelease-1>', set_eps12)
    slider_eps12.set(target_eps12)
    slider_eps12.grid(row=2,column=1)
    Label(frame, textvariable=value_eps12).grid(row=2, column=2)
    frame.pack()

    tkroot.update()

# wrapper on GnuPlot via Pizza.py gnu tool

if me == 0:
    from gnu import gnu

    gn = gnu()
    gn.plot(xaxis, yaxis)
    gn.title(compute, "Timestep", "Temperature")

# endless loop, checking status of GUI settings every Nfreq steps
# run with pre yes/no and post yes/no depending on go/stop status
# re-invoke fix langevin with new seed when temperature slider changes
# after re-invoke of fix langevin, run with pre yes

running = 0
# define current values
temp = target_temp
eps11 = target_eps11
eps12 = target_eps12
eps22 = target_eps22

seed = 12345

lmp.command("fix temp all temp/csvr %g %g 0.05 %d" % (temp, temp, seed))
# lmp.command("fix temp all temp/berendsen %g %g 0.05 %d" % (temp,temp,seed))
# lmp.command("fix md all nvt temp %g %g 0.05" % (temp,temp,seed))

while 1:
    if me == 0: tkroot.update()

    # communication of calues (from root)
    running, runflag, breakflag = comm.bcast([running, runflag, breakflag], root=0)
    temp, target_temp, seed = comm.bcast([temp, target_temp, seed], root=0)
    eps11, target_eps11 = comm.bcast([eps11, target_eps11], root=0)
    eps12, target_eps12 = comm.bcast([eps12, target_eps12], root=0)
    eps22, target_eps22 = comm.bcast([eps22, target_eps22], root=0)

    # handle sliders
    if temp != target_temp:
        temp = target_temp
        seed += me + 1
        lmp.command("fix temp all temp/csvr %g %g 0.05 %d" % (temp, temp, seed))
        # lmp.command("fix temp all temp/berendsen %g %g 0.05 %d" % (temp,temp,seed))
        # lmp.command("fix md all nvt temp %g %g 0.05" % (temp,temp,seed))
        running = 0
    if eps11 != target_eps11:
        lmp.command(str_pair_coeff_11 % target_eps11)
        running = 0
    if eps12 != target_eps12:
        lmp.command(str_pair_coeff_12 % target_eps12)
        running = 0
    if eps22 != target_eps22:
        lmp.command(str_pair_coeff_22 % target_eps22)
        running = 0

    # handle runflags
    if runflag and running:
        lmp.command("run %d pre no post no" % nfreq)
        ntimestep += nfreq
        if me == 0: update(ntimestep)
    elif runflag and not running:
        lmp.command("run %d pre yes post no" % nfreq)
        ntimestep += nfreq
        if me == 0: update(ntimestep)
    elif not runflag and running:
        lmp.command("run %d pre no post yes" % nfreq)
        ntimestep += nfreq
        if me == 0: update(ntimestep)
    if breakflag:
        break
    if runflag:
        running = 1
    else:
        running = 0
    time.sleep(0.01)

lmp.command("run 0 pre no post yes")

print("Proc %d out of %d procs has" % (me, nprocs) + str(lmp))
lmp.close()
