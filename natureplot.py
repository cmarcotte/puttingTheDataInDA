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

times = range(1,1001)

stateDir = f"{workDir}/state/"

hx = 0.015

print("Writing images:")
widgets = [progressbar.Percentage(), progressbar.ETA(), progressbar.Bar()]
bar = progressbar.ProgressBar(widgets=widgets, max_value=len(times)).start()

fig, axs = plt.subplots(3, 3, figsize=(sw,sw), sharex="col", sharey="row", constrained_layout=True)

# This figure is a 3x3 grid, showing slices through U,V,W (columns) across x,y,z (rows)

X = 200
Y = 200
Z =  50
x = 100 #np.random.randint(0,high=X)
y = 100 #np.random.randint(0,high=Y)
z = 25 #np.random.randint(0,high=Z)

cmaps = ["Oranges", "Purples", "Greens"]
varnames=["u", "v", "w"]

for n,t in enumerate(times):
	try:
		S 	= readData.readRestarts([f"{stateDir}/{t:04d}/restart3d.{n:03d}" for n in range(4)])
		for o in range(3):
			for m in range(3):
				axs[o,m].cla()
				axs[o,m].axes.set_aspect('equal')
			axs[o,0].pcolormesh(np.arange(0,50)*hx, np.arange(0,200)*hx, S[o][x,:,:], snap=True, shading="auto", rasterized=True, vmin=0.0, vmax=1.0, cmap=cmaps[o])
			axs[o,1].pcolormesh(np.arange(0,50)*hx, np.arange(0,200)*hx, S[o][:,y,:], snap=True, shading="auto", rasterized=True, vmin=0.0, vmax=1.0, cmap=cmaps[o])
			axs[o,2].pcolormesh(np.arange(0,200)*hx, np.arange(0,200)*hx, S[o][:,:,z], snap=True, shading="auto", rasterized=True, vmin=0.0, vmax=1.0, cmap=cmaps[o])
			if n==0:
				im = axs[o,2].pcolormesh(np.arange(0,200)*hx, np.arange(0,200)*hx, S[o][:,:,z], snap=True, shading="auto", rasterized=True, vmin=0.0, vmax=1.0, cmap=cmaps[o])
				clb = plt.colorbar(im, ax=axs[o,:])
				clb.ax.set_title(r'${0}$'.format(varnames[o]))
			for m in range(3):
				axs[o,m].set_xticks([])
				axs[o,m].set_yticks([])
	#	axs[0,0].set_title(r'$x={0}$ [cm]'.format(np.round(x*hx,decimals=2)))
	#	axs[0,1].set_title(r'$y={0}$ [cm]'.format(np.round(y*hx,decimals=2)))
	#	axs[0,2].set_title(r'$z={0}$ [cm]'.format(np.round(z*hx,decimals=2)))
		axs[0,0].set_title(r'$x={0}/{1}$'.format(x,X))
		axs[0,1].set_title(r'$y={0}/{1}$'.format(y,Y))
		axs[0,2].set_title(r'$z={0}/{1}$'.format(z,Z))
		plt.suptitle(r"$t={0}$ [ms]".format(2.0*t))
		plt.savefig(f"./continuity_frames/{t:04d}.png", dpi=300, bbox_inches="tight")
		bar.update(n + 1)
	except:
		pass
bar.finish()
