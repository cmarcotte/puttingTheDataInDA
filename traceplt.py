import numpy as np
import progressbar
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc

plt.style.use('seaborn-paper')

rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)
matplotlib.rcParams['axes.titlesize'] = 10
matplotlib.rcParams['axes.labelsize'] = 10
matplotlib.rcParams['xtick.labelsize'] = 9
matplotlib.rcParams['ytick.labelsize'] = 9

sw = 3.40457
dw = 7.05826

import argparse

parser = argparse.ArgumentParser(description="Load and visualize restart / obs files for a given input directory.")
parser.add_argument('--workDir', type=str, help="Data base directory.")
args = parser.parse_args()
workDir = f"/data/{args.workDir}"

import readData

guesDir = f"{workDir}/gues/mean/"
analDir = f"{workDir}/anal/mean/"
obsDir  = f"{workDir}/obs/"

import os
n_dirs = 0
for root, dirs, files in os.walk(guesDir, topdown=True):
        n_dirs += len(dirs)
times = range(1,n_dirs+1)
times = range(1,1001)

# pre-allocate (are these going to be written to or just over-written and re-allocated? python...)
gU = np.zeros((200,200,50))
gV = np.zeros((200,200,50))
gW = np.zeros((200,200,50))
aU = np.zeros((200,200,50))
aV = np.zeros((200,200,50))
aW = np.zeros((200,200,50))
oU = np.zeros((200,200,50))

hx = 0.015

def obs2arrs(Obs):
        tmp = np.vstack(Obs)
        x = tmp[:,3]-1
        y = tmp[:,2]-1
        z = tmp[:,1]-1
        u = tmp[:,4]
        return x,y,z,u


print("Computing Errors from Observation and Background Ensemble:")
widgets = [progressbar.Percentage(), progressbar.ETA(), progressbar.Bar()]
bar = progressbar.ProgressBar(widgets=widgets, max_value=len(times)).start()

trcs = np.NaN*np.zeros((len(times), 7))
for n,t in enumerate(times[0:1]):
	Obs		= readData.readObservations(f"{obsDir}/{t:04d}.dat")
	X,Y,Z,U 	= obs2arrs(Obs)
	X = X.astype(int); Y = Y.astype(int); Z = Z.astype(int);
	oU *= np.NaN
	oU[X,Y,Z] = U

oind = np.argmin((X-99)**2 + (Y-159)**2 + (Z-0)**2) #np.random.randint(0,high=len(U))

x = X[oind]
y = Y[oind]
z = Z[oind]

for n,t in enumerate(times):
	try:
		gU,gV,gW	= readData.readRestarts([f"{guesDir}/{t:04d}/restart3d.{n:03d}" for n in range(4)])
		aU,aV,aW	= readData.readRestarts([f"{analDir}/{t:04d}/restart3d.{n:03d}" for n in range(4)])
		Obs		= readData.readObservations(f"{obsDir}/{t:04d}.dat")
		X,Y,Z,U 	= obs2arrs(Obs)
		X = X.astype(int); Y = Y.astype(int); Z = Z.astype(int);
		oU *= np.NaN
		oU[X,Y,Z] = U
		trcs[n,0] = gU[X[oind],Y[oind],Z[oind]]
		trcs[n,1] = gV[X[oind],Y[oind],Z[oind]]
		trcs[n,2] = gW[X[oind],Y[oind],Z[oind]]
		trcs[n,3] = aU[X[oind],Y[oind],Z[oind]]
		trcs[n,4] = aV[X[oind],Y[oind],Z[oind]]
		trcs[n,5] = aW[X[oind],Y[oind],Z[oind]]
		trcs[n,6] = U[oind]
	except:
		pass
	bar.update(n + 1)
bar.finish()

fig, axs = plt.subplots(3, 1, figsize=(sw, sw), sharex=True, sharey=True, constrained_layout=True)
labs = ["(a)", "(b)", "(c)"]
vrns = ["$u$", "$v$", "$w$"]
t = 2.0*np.array(times)
for n in range(3):
	axs[n].plot(t, trcs[:,n], f"-C0", label=r"{0}({1},{2},{3})".format("Background ",x,y,z))
	axs[n].plot(t, trcs[:,3+n], f"-C1", label=r"{0}({1},{2},{3})".format("Analysis ",x,y,z))
	if n==0:
		axs[n].plot(t, trcs[:,-1], ".-C2", linewidth=0.5, markersize=4, label="Obs.")
	axs[n].set_ylabel(r"{0}".format(vrns[n]))
	axs[n].text(-0.2, +1.1, labs[n], horizontalalignment="center", verticalalignment="center", transform=axs[n].transAxes)
axs[0].legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False)
axs[-1].set_xlabel(r"$t$ [ms]")
axs[-1].set_xlim([0,2000])
#axs[0].set_ylim([0.0,1.0])
#axs[0].set_yscale("log")
plt.savefig(f"{args.workDir}_traces.pdf", bbox_inches="tight", dpi=600)
plt.close()
