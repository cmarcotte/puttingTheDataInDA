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

meanDir = f"{workDir}/gues/mean/"
sprdDir = f"{workDir}/gues/sprd/"
obsDir  = f"{workDir}/obs/"

import os
n_dirs = 0
for root, dirs, files in os.walk(meanDir, topdown=True):
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

print("Writing images:")
widgets = [progressbar.Percentage(), progressbar.ETA(), progressbar.Bar()]
bar = progressbar.ProgressBar(widgets=widgets, max_value=len(times)).start()

fig, axs = plt.subplots(3, 3, figsize=(sw,sw), sharex=True, sharey=True, constrained_layout=True)

for n,t in enumerate(times):
	for o in range(3):
		for m in range(3):
			axs[o,m].cla()
			axs[o,m].axes.set_aspect('equal')
	try:
		oU *= np.NaN
		gU,_,_ 	= readData.readRestarts([f"{meanDir}/{t:04d}/restart3d.{n:03d}" for n in range(4)])
		sU,_,_ 	= readData.readRestarts([f"{sprdDir}/{t:04d}/restart3d.{n:03d}" for n in range(4)])
		Obs	= readData.readObservations(f"{obsDir}/{t:04d}.dat")
		#for m,obs in enumerate(Obs):
		#	oU[int(obs[3])-1,int(obs[2])-1,int(obs[1])-1] = obs[4]
		X,Y,Z,U = obs2arrs(Obs)
		for o,z in enumerate([0,24,49]):
			axs[o,0].pcolormesh(np.arange(0,200)*hx, np.arange(0,200)*hx, gU[:,:,z], snap=True, shading="auto", rasterized=True, vmin=0.0, vmax=1.0, cmap="Oranges")
			axs[o,1].pcolormesh(np.arange(0,200)*hx, np.arange(0,200)*hx, sU[:,:,z], snap=True, shading="auto", rasterized=True, vmin=0.0, vmax=1.0, cmap="Oranges")
			#axs[o,2].pcolormesh(np.arange(0,200)*hx, np.arange(0,200)*hx, oU[:,:,z], snap=True, shading="auto", rasterized=True, vmin=0.0, vmax=1.0, cmap="Oranges")
			iz = np.flatnonzero(Z==z)
			if len(iz) > 0:
				# should these be (X+1,Y+1) because the array formation subtracts one to make (X,Y) work as 0-based indices?
				axs[o,2].tricontourf((1+Y[iz])*hx, (1+X[iz])*hx, U[iz], levels=129, vmin=0.0, vmax=1.0, cmap="Oranges")
				if len(iz) < 100:
					axs[o,2].plot((1+Y[iz])*hx, (1+X[iz])*hx, ".k", markersize=3)
			#axs[o,2].contourf(np.arange(0,200)*hx, np.arange(0,200)*hx, oU[:,:,z], vmin=0.0, vmax=1.0, cmap="Oranges")
			#axs[o,0].imshow(gU[:,:,z], 		interpolation="none", vmin=0.0, vmax=1.0, cmap="Oranges")
			#axs[o,1].imshow(aU[:,:,z], 		interpolation="none", vmin=0.0, vmax=1.0, cmap="Oranges")
			#axs[o,2].imshow(oU[0:-1:3,0:-1:3,z], 	interpolation="none", vmin=0.0, vmax=1.0, cmap="Oranges")
			if z==0:
				obsdepth = 0.0
			else:
				obsdepth = (z+1)*hx
			axs[o,0].set_ylabel(r"${0}$ [cm]".format(obsdepth))
			for m in range(3):
				#axs[o,m].axis("off")
				axs[o,m].set_xticks([])
				axs[o,m].set_yticks([])
		if n==0:
			im = axs[o,0].pcolormesh(np.arange(0,200)*hx, np.arange(0,200)*hx, gU[:,:,0], snap=True, shading="auto", rasterized=True, vmin=0.0, vmax=1.0, cmap="Oranges")
			clb = plt.colorbar(im, ax=axs)
			clb.ax.set_title(r'$u$')
			axs[0,0].set_xlim([0*hx,199*hx])
			axs[0,0].set_ylim([1*hx,200*hx])
		axs[0,0].set_title(f"Mean")
		axs[0,1].set_title(f"Spread")
		axs[0,2].set_title(f"Obs.")
		plt.suptitle(r"$t={0}$ [ms]".format(2.0*t))
		plt.savefig(f"./dynamics_frames/{t:04d}.png", dpi=300, bbox_inches="tight")
		bar.update(n + 1)
	except:
		pass
bar.finish()
