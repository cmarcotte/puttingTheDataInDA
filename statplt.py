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
sprdDir = f"{workDir}/gues/sprd/"
obsDir  = f"{workDir}/obs/"

import os
n_dirs = 0
for root, dirs, files in os.walk(guesDir, topdown=True):
        n_dirs += len(dirs)
times = range(1,n_dirs+1)
times = range(1,1001)

# pre-allocate (are these going to be written to or just over-written and re-allocated? python...)
gU = np.zeros((200,200,50))
sU = np.zeros((200,200,50))
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

errs = np.NaN*np.zeros((len(times), 200, 200, 50,))
sprd = np.NaN*np.zeros((len(times), 200, 200, 50,))
for n,t in enumerate(times):
	try:
		gU,_,_	= readData.readRestarts([f"{guesDir}/{t:04d}/restart3d.{n:03d}" for n in range(4)])
		sU,_,_	= readData.readRestarts([f"{sprdDir}/{t:04d}/restart3d.{n:03d}" for n in range(4)])
		Obs	= readData.readObservations(f"{obsDir}/{t:04d}.dat")
		X,Y,Z,U = obs2arrs(Obs)
		X = X.astype(int); Y = Y.astype(int); Z = Z.astype(int);
		oU *= np.NaN
		oU[X,Y,Z] = U
		errs[n,X,Y,Z] = np.abs(gU[X,Y,Z]-oU[X,Y,Z])
		sprd[n,X,Y,Z] = sU[X,Y,Z]
	except:
		pass
	bar.update(n + 1)
bar.finish()

fig, axs = plt.subplots(2, 1, figsize=(sw, sw), sharex=True, sharey=True, constrained_layout=True)
labels = ["EPI", "ENDO"]
labs = ["(a)", "(b)"]
Zs = [0,49]; #np.unique(Z)
for o,z in enumerate(Zs):
	x = 2.0*np.array(times)
	iz = np.flatnonzero(Z==z)
	y = np.sqrt(np.mean(errs[:,X[iz],Y[iz],Z[iz]]**2,axis=-1))
	dy= np.std(errs[:,X[iz],Y[iz],Z[iz]],axis=-1)
	sy= np.sqrt(np.mean(sprd[:,X[iz],Y[iz],Z[iz]]**2,axis=-1))
	ds = np.std(sprd[:,X[iz],Y[iz],Z[iz]],axis=-1)
	axs[o].fill_between(x,  y - dy,  y + dy, color="C0", alpha=0.3, label="")
	axs[o].fill_between(x, sy - ds, sy + ds, color="C1", alpha=0.3, label="")
	axs[o].plot(x, y, f"-C0", label="RMS Error")
	axs[o].plot(x, sy, "-C1", linewidth=1, label="RMS Spread")
	axs[o].set_ylabel(r"{0}".format(labels[o]))
	axs[o].text(-0.2, +1.1, labs[o], horizontalalignment="center", verticalalignment="center", transform=axs[o].transAxes)
axs[0].legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=False)
axs[-1].set_xlabel(r"$t$ [ms]")
axs[-1].set_xlim([0,2000])
axs[0].set_ylim([0.0,1.0])
#axs[0].set_yscale("log")
plt.savefig(f"{args.workDir}_SurfaceErrors.pdf", bbox_inches="tight", dpi=600)
plt.close()
