import numpy as np
import progressbar
import matplotlib.pyplot as plt

# these are consts which may need to be adapted to your problem
nx,ny,nz,nnz = 50,200,200,50

def readRestarts(fnames, nproc=4):
	"""
	Read (u,v,w) fields from Fortran data files
		fnames: list of restart filenames (including path)
	"""
	u = np.zeros((nz+0,ny+0,nx+0))
	v = np.zeros((nz+0,ny+0,nx+0))
	w = np.zeros((nz+0,ny+0,nx+0))
	for p in range(nproc):
		with open(fnames[p],'rb') as fid:
			garbage = np.fromfile(fid, np.float32, 1)
			u_tmp = np.fromfile(fid, np.float64, (nx+2)*(ny+2)*(nnz+2))
			v_tmp = np.fromfile(fid, np.float64, (nx+1)*(ny+1)*(nnz+1))
			w_tmp = np.fromfile(fid, np.float64, (nx+1)*(ny+1)*(nnz+1))
			u_tmp = np.reshape(u_tmp,(nnz+2,ny+2,nx+2))
			u_tmp = u_tmp[1:-1,1:-1,1:-1]
			v_tmp = np.reshape(v_tmp,(nnz+1,ny+1,nx+1))
			v_tmp = v_tmp[1:,1:,1:]
			w_tmp = np.reshape(w_tmp,(nnz+1,ny+1,nx+1))
			w_tmp = w_tmp[1:,1:,1:]
			u[p*nnz:(p+1)*nnz,:,:] = u_tmp
			v[p*nnz:(p+1)*nnz,:,:] = v_tmp
			w[p*nnz:(p+1)*nnz,:,:] = w_tmp
	return (u,v,w)

def readTips(fnames, nproc=4):
	"""
	Read tip files from simulation output
		fnames: list of tip filenames (including path)
	"""
	for p in range(nproc):
		tmp = np.loadtxt(fnames[p])
		if tmp.size != 0:
			try:
				tipz = np.concatenate((tipz, tmp), axis=0)
			except:
				tipz = tmp
	return tipz

def readObservations(fname):
	"""
	Read observation files generated with surfobsmake.sh
		fnames: observation filename (including path) NOTE: not a list!
	"""
	with open(fname,'rb') as fid:
		EOF = False
		while not EOF:
			garbage = np.fromfile(fid, np.float32, 1)
			tmp = np.fromfile(fid, np.float32, 6)
			garbage = np.fromfile(fid, np.float32, 1)
			if len(tmp)==0:
				EOF = True
			else:
				EOF = False
				tmp = np.reshape(np.array(tmp),(1,6))
				try:
					obs = np.concatenate((obs, tmp), axis=0)
				except:
					obs = tmp
	return obs

def collectStates(basedir, times, plotting=False):
	U = np.zeros((len(times), nz, ny, nx))
	V = np.zeros((len(times), nz, ny, nx))
	W = np.zeros((len(times), nz, ny, nx))
	print("Collecting states:")
	widgets = [progressbar.Percentage(), progressbar.Bar()]
	bar = progressbar.ProgressBar(widgets=widgets, max_value=len(times)).start()
	for (nt,t) in enumerate(times):
		# state file names
		fnames = [f'{basedir}/{t:04d}/restart3d.{n:03d}' for n in range(4)]
		u,v,w = readRestarts(fnames)
		# write into collection arrays
		U[nt] = u
		V[nt] = v
		W[nt] = w
		if plotting:
			fig, axs = plt.subplots(1,3,sharex=True, sharey=True, constrained_layout=True)
			axs[0].imshow(u[:,:,25], vmin=0.0, vmax=1.0)
			axs[1].imshow(v[:,:,25], vmin=0.0, vmax=1.0)
			axs[2].imshow(w[:,:,25], vmin=0.0, vmax=1.0)
			plt.savefig(f"./{t:04d}.png", bbox_inches='tight')
			plt.close()
		bar.update(nt + 1)
	bar.finish()
	return U,V,W

def collectObservations(basedir, times):
	O = []
	print("Collecting Observations:")
	widgets = [progressbar.Percentage(), progressbar.Bar()]
	bar = progressbar.ProgressBar(widgets=widgets, max_value=len(times)).start()
	for (nt,t) in enumerate(times):
		# observation file name
		fname = f'{basedir}/{t:04d}.dat'
		obs = readObservations(fname)
		# write into collection arrays
		O.append(obs)
		bar.update(nt + 1)
	bar.finish()
	return O

def collectTips(basedir, times):
	T = []
	print("Collecting Tips:")
	widgets = [progressbar.Percentage(), progressbar.Bar()]
	bar = progressbar.ProgressBar(widgets=widgets, max_value=len(times)).start()
	for (nt,t) in enumerate(times):
		# tip file names
		fnames = [f'{basedir}/{t:04d}/tips{n:03d}.0001' for n in range(4)]
		tips = readTips(fnames)
		# write into collection arrays
		T.append(tips)
		bar.update(nt + 1)
	bar.finish()
	return T

if __name__ == "__main__":
	# example data; your usage will vary
	stateDir	= "/home/chris/Development/Fortran/cardiacDA/3D/DATA/nature_4cm_straight_60degree_2/"
	stateTimes 	= range(0,10)
	obsDir		= "/home/chris/Development/Alessio Data/obs/"
	obsTimes 	= range(140,150)
	tipsDir	= "/home/chris/Development/Fortran/cardiacDA/3D/DATA/nature_2D/tips/"
	tipsTimes	= range(141,150)
	# example usage for single read
	u,v,w 		= readRestarts([f'{stateDir}/{stateTimes[0]:04d}/restart3d.{n:03d}' for n in range(4)])
	#print((u,v,w))
	obs 		= readObservations(f'{obsDir}/{obsTimes[0]:04d}.dat')
	print(obs)
	tips		= readTips([f'{tipsDir}/{tipsTimes[0]:04d}/tips{n:03d}.0001' for n in range(4)])
	print(tips)
	# example usage for sequence read
	U,V,W 		= collectStates(stateDir, stateTimes, plotting=False)	# plotting=True is much slower
	Obs		= collectObservations(obsDir, obsTimes)
	Tips		= collectTips(tipsDir, tipsTimes)
