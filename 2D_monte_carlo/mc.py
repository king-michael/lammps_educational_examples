#!/usr/bin/env python -i
# preceeding line should have path for Python on your machine

# mc.py
# Purpose: mimic operation of example/MC/in.mc via Python
# Syntax:  mc.py in.mc
#          in.mc = LAMMPS input script

import sys,random,math,time
sys.path.append("../pizza")
# set these parameters
# make sure neigh skin (in in.mc) > 2*deltamove

#nloop = 3000
deltaperturb = 0.2
deltamove = 0.1
kT = 0.05
random.seed(27848)



# methods called by GUI

def run():
  global runflag
  runflag = 1
def stop():
  global runflag
  runflag = 0
def restart():
  global runflag,x
  runflag = 0
  lmp.command("undump python")
  lmp.command("dump python all atom %d tmp.dump" % nfreq)
  lmp.command("run 0")
  x = lmp.extract_atom("x",3)
  #initial()
  unorder()
def setdeltamove(value):
  global deltamovetarget
  deltamovetarget = slider.get()
def quit():
  global breakflag
  breakflag = 1


# method called by timestep loop every Nfreq steps
# read dump snapshot and viz it, update plot with compute value

def update(ntimestep):
  d.next()
  d.unscale()
  p.single(ntimestep)
  pm.load("tmp.pdb")
  pm.forward()
  value = lmp.extract_compute("thermo_pe",0,0)
  xaxis.append(ntimestep)
  yaxis.append(value)
  gn.plot(xaxis,yaxis)

# parse command line

argv = sys.argv
if len(argv) != 2:
  print "Syntax: vizplotgui_pymol.py in.lammps Nfreq compute-ID"
  sys.exit()

infile = sys.argv[1]
nfreq =  1  # int(sys.argv[2])
#compute = sys.argv[3]

me = 0

# uncomment if running in parallel via Pypar
#import pypar
#me = pypar.rank()
#nprocs = pypar.size()




from lammps import lammps
lmp = lammps()

# run infile one line at a time
# just sets up MC problem

# run 0 to get energy of perfect lattice
# emin = minimum energy


def initial():
	global x, lines, ntimestep, deltamovetarget
	global value, xaxis, yaxis
	global breakflag,runflag
	global e, elast, kT, natoms, estart, emin
	global deltamove, naccept, deltaperturb
	lmp.command("clear")
	lines = open(infile,'r').readlines()
	for line in lines: lmp.command(line)
	lmp.command("variable e equal pe")

	lmp.command("thermo %d" % nfreq)
	lmp.command("dump python all atom %d tmp.dump" % nfreq)
	lmp.command("run 0")
	value = lmp.extract_compute("thermo_pe",0,0)
	ntimestep = 0
	xaxis = [ntimestep]
	yaxis = [value]

	breakflag = 0
	runflag = 0
	deltamovetarget = deltamove


	natoms = lmp.extract_global("natoms",0)
	emin = lmp.extract_compute("thermo_pe",0,0) / natoms
	lmp.command("variable emin equal $e")
	x = lmp.extract_atom("x",3)
initial()
# disorder the system
# estart = initial energy


def unorder():
	global x, lines, ntimestep, deltamovetarget
	global value, xaxis, yaxis
	global breakflag,runflag
	global e, elast, kT, natoms, estart, emin
	global deltamove, naccept, deltaperturb
	for i in xrange(natoms):
		x[i][0] += deltaperturb * (2*random.random()-1)
		x[i][1] += deltaperturb * (2*random.random()-1)

	lmp.command("variable elast equal $e")
	lmp.command("thermo_style custom step v_emin v_elast pe")
	lmp.command("run 0")
	x = lmp.extract_atom("x",3)
	lmp.command("variable elast equal $e")
	
	estart = lmp.extract_compute("thermo_pe",0,0) / natoms

	# loop over Monte Carlo moves
	# extract x after every run, in case reneighboring changed ptr in LAMMPS

	elast = estart
	naccept = 0
unorder()

# wrapper on PyMol
# just proc 0 handles reading of dump file and viz

if me == 0:
  import pymol
  pymol.finish_launching()

  from dump import dump
  from pdbfile import pdbfile
  from pymol import cmd as pm

  d = dump("tmp.dump",0)
  p = pdbfile(d)
  d.next()
  d.unscale()
  p.single(ntimestep)
  pm.load("tmp.pdb")
  pm.show("spheres","tmp")
  #pm.zoom()
  pm.set_view("(\
     1.000000000,    0.000000000,    0.000000000,\
     0.000000000,    1.000000000,    0.000000000,\
     0.000000000,    0.000000000,    1.000000000,\
     0.000000000,    0.000000000,  -30.770051956,\
     5.363590240,    4.751065731,    0.000000000,\
    24.259342194,   37.280761719,  -20.000000000 )")
# display GUI with run/stop buttons and slider for temperature

if me == 0:
  from Tkinter import *
  tkroot = Tk()
  tkroot.withdraw()
  root = Toplevel(tkroot)
  root.title("LAMMPS GUI")

  frame = Frame(root)
  Button(frame,text="Run",command=run).pack(side=LEFT)
  Button(frame,text="Stop",command=stop).pack(side=LEFT)
  #Button(frame,text="Restart",command=restart).pack(side=LEFT)
  slider = Scale(frame,from_=0.01,to=0.5,resolution=0.01,
                 orient=HORIZONTAL,label="deltamove")
  slider.bind('<ButtonRelease-1>',setdeltamove)
  slider.set(deltamovetarget)
  slider.pack(side=LEFT)
  Button(frame,text="Quit",command=quit).pack(side=RIGHT)
  frame.pack()
  tkroot.update()

# wrapper on GnuPlot via Pizza.py gnu tool

if me == 0:
  from gnu import gnu
  gn = gnu()
  gn.plot(xaxis,yaxis)
  gn.title("thermo_pe","Timestep","PotEng")

# endless loop, checking status of GUI settings every Nfreq steps
# run with pre yes/no and post yes/no depending on go/stop status
# re-invoke fix langevin with new seed when temperature slider changes
# after re-invoke of fix langevin, run with pre yes

running = 0
deltamove = deltamovetarget
seed = 12345

def domc():
  global x, lines, ntimestep, deltamovetarget
  global value, xaxis, yaxis
  global breakflag,runflag
  global e, elast, kT, natoms, estart, emin
  global deltamove, naccept, deltaperturb
  iatom = random.randrange(0,natoms)
  x0 = x[iatom][0]
  y0 = x[iatom][1]

  x[iatom][0] += deltamove * (2*random.random()-1)
  x[iatom][1] += deltamove * (2*random.random()-1)

  lmp.command("run 1 pre no post no")
  x = lmp.extract_atom("x",3)
  e = lmp.extract_compute("thermo_pe",0,0) / natoms

  if e <= elast:
    elast = e
    lmp.command("variable elast equal $e")
    naccept += 1
  elif random.random() <= math.exp(natoms*(elast-e)/kT):
    elast = e
    lmp.command("variable elast equal $e")
    naccept += 1
  else:
    x[iatom][0] = x0
    x[iatom][1] = y0
    
    
    
while 1:
  if me == 0: tkroot.update()
  if deltamove != deltamovetarget:
    deltamove = deltamovetarget
    seed += me+1
    running = 0
  if runflag and running:
    domc()
    ntimestep += nfreq
    if me == 0: update(ntimestep)
  elif runflag and not running:
    domc()
    ntimestep += nfreq
    if me == 0: update(ntimestep)
  elif not runflag and running:
    domc()
    ntimestep += nfreq
    if me == 0: update(ntimestep)
  if breakflag: break
  if runflag: running = 1
  else: running = 0
  time.sleep(0.01)


    
# final energy and stats

lmp.command("variable nbuild equal nbuild")
nbuild = lmp.extract_variable("nbuild",None,0)

lmp.command("run 0")
estop = lmp.extract_compute("thermo_pe",0,0) / natoms

print "MC stats:"
print "  starting energy =",estart
print "  final energy =",estop
print "  minimum energy of perfect lattice =",emin
print "  accepted MC moves =",naccept
print "  neighbor list rebuilds =",nbuild


