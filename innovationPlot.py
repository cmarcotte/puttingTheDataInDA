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

import os
n_dirs = 0
for root, dirs, files in os.walk(guesDir, topdown=True):
	n_dirs += len(dirs)
times = range(1,n_dirs+1)
times = range(1,1001)

# pre-allocate (are these going to be written to or just over-written and re-allocated? python...)
gU = np.zeros((200,200,50))
aU = np.zeros((200,200,50))

hx = 0.015

print("Writing images:")
widgets = [progressbar.Percentage(), progressbar.ETA(), progressbar.Bar()]
bar = progressbar.ProgressBar(widgets=widgets, max_value=len(times)).start()

fig, axs = plt.subplots(3, 6, figsize=(dw,sw), sharex=True, sharey=True, constrained_layout=True)

for n,t in enumerate(times):
	for o in range(6):
		for m in range(3):
			axs[m,o].cla()
			axs[m,o].axes.set_aspect('equal')
	try:
		gU,_,_ 	= readData.readRestarts([f"{guesDir}/{t:04d}/restart3d.{q:03d}" for q in range(4)])
		aU,_,_ 	= readData.readRestarts([f"{analDir}/{t:04d}/restart3d.{q:03d}" for q in range(4)])
		for o,z in enumerate([0,10,20,30,40,49]):
			axs[0,o].pcolormesh(np.arange(0,200)*hx, np.arange(0,200)*hx, gU[:,:,z], 		snap=True, shading="auto", rasterized=True, vmin=0.0, vmax=1.0, cmap="Oranges")
			axs[1,o].pcolormesh(np.arange(0,200)*hx, np.arange(0,200)*hx, aU[:,:,z]-gU[:,:,z], 	snap=True, shading="auto", rasterized=True, vmin=-0.1, vmax=+0.1, cmap="seismic")
			axs[2,o].pcolormesh(np.arange(0,200)*hx, np.arange(0,200)*hx, aU[:,:,z], 		snap=True, shading="auto", rasterized=True, vmin=0.0, vmax=1.0, cmap="Oranges")
			if z==0:
				obsdepth = 0.0
			else:
				obsdepth = (z+1)*hx
			axs[-1,o].set_xlabel(r"${0:2.3f}$ [cm]".format(obsdepth))
			for m in range(3):
				axs[m,o].set_xticks([])
				axs[m,o].set_yticks([])
		if n==0:
			im = axs[0,o].pcolormesh(np.arange(0,200)*hx, np.arange(0,200)*hx, gU[:,:,z], 		snap=True, shading="auto", rasterized=True, vmin=0.0, vmax=1.0, cmap="Oranges")
			clb = plt.colorbar(im, ax=axs[0,:])
			clb.ax.set_ylabel(r'$u^b$')
			im = axs[1,o].pcolormesh(np.arange(0,200)*hx, np.arange(0,200)*hx, aU[:,:,z]-gU[:,:,z],	snap=True, shading="auto", rasterized=True, vmin=-0.1, vmax=+0.1, cmap="seismic")
			clb = plt.colorbar(im, ax=axs[1,:])
			clb.ax.set_ylabel(r'$u^a-u^b$')
			im = axs[2,o].pcolormesh(np.arange(0,200)*hx, np.arange(0,200)*hx, gU[:,:,z], 		snap=True, shading="auto", rasterized=True, vmin=0.0, vmax=1.0, cmap="Oranges")
			clb = plt.colorbar(im, ax=axs[2,:])
			clb.ax.set_ylabel(r'$u^a$')
			axs[0,0].set_xlim([0*hx,199*hx])
			axs[0,0].set_ylim([1*hx,200*hx])
		plt.suptitle(r"$t={0}$ [ms]".format(2.0*t))
		plt.savefig(f"./innovation_frames/{t:04d}.png", dpi=300, bbox_inches="tight")
		bar.update(n + 1)
	except:
		pass
bar.finish()
