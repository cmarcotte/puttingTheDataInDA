include("loadAlessioData.jl")
using .loadAlessioData, DelimitedFiles, FortranFiles, Printf, ProgressBars, PyPlot, Peaks, Dierckx, Interpolations

function main(; plotting=false, movie=false, outputObs=true, l=3.0)
	
	# where is the data
	basepath = "./2011-05-28_Rec78-103_Pace_Apex/"
	inds = getPossibleIndices(basepath)
	ind = 114 #rand(inds) # This is just useful for debugging; normally you would choose an index
	
	# print some diagnostic notices
	print("Loading data for $(ind)...\t")
	
	# load in some data 
	states = loadFile(basepath, ind)
	
	print("Done.\n")
	
	# form the space and time masks
	spaceMasks = [spatialMask(state) for state in states]
	timeMasks  = [temporalMask(state) for state in states]
	mutualSpaceMask = mutualMask(spaceMasks)
	mutualTimeMask = mutualMask(timeMasks)
	
	# determine mutualSpaceMask centroid for square center
	center = spatialMaskCentroid(mutualSpaceMask)
	corners = loadAlessioData.square(theta=-0.0*pi/18, l=l, center=center)
	
	# make some plots
	if plotting
	# sampling position is randomly sampled and checked against mutual space mask
		samplePoint = Int.(center./0.06) #[rand(1:128), rand(1:128)]
		while mutualSpaceMask[samplePoint[1],samplePoint[2]] == false
			samplePoint = [rand(1:128), rand(1:128)]
		end
		t, V = sampleData(states, samplePoint, mutualTimeMask, mutualSpaceMask)
		fig, axs = plt.subplots(3,1,figsize=(loadAlessioData.sw,loadAlessioData.dw),constrained_layout=true)
		plotTrace!(fig, axs, t, V; threshold=0.25)
		plt.savefig("./trace/$(ind)_trace.svg", bbox_inches="tight")
		plt.savefig("./trace/$(ind)_trace.pdf", bbox_inches="tight")
		plt.close()
	
		# plotting a single state at a prescribed time
		sampleInt = argmax(abs.(V[1].-V[2]))
		#sampleTime = rand((1:length(mutualTimeMask))[mutualTimeMask])
		fig, axs = plt.subplots(2,1,figsize=(loadAlessioData.sw,loadAlessioData.dw),constrained_layout=true)
		plotState!(fig, axs, sampleInt, states, spaceMasks, timeMasks, [size(state) for state in states], mutualSpaceMask, mutualTimeMask; drawColorbar=true)
		plotSamplingSquare!(fig, axs, corners)
		plt.savefig("./state/$(ind)_state.svg", bbox_inches="tight")
		plt.savefig("./state/$(ind)_state.pdf", bbox_inches="tight")
		plt.close()
		
		@show sampleInt
		@show t[sampleInt]
		
		# making a figure for the paper
		fig = plt.figure(figsize=(loadAlessioData.sw,loadAlessioData.sw), constrained_layout=true)
		gs = plt.GridSpec(2, 2, figure=fig)
		ax1 = plt.subplot(gs.new_subplotspec((1, 0), colspan=2))
		ax4 = plt.subplot(gs.new_subplotspec((0, 0)))
		ax5 = plt.subplot(gs.new_subplotspec((0, 1)))
		axs = [ax4,ax5,ax1]
		plotState!(fig, axs, sampleInt, states, spaceMasks, timeMasks, [size(state) for state in states], mutualSpaceMask, mutualTimeMask; drawColorbar=true)
		plotSamplingSquare!(fig, axs, corners)
		axs[1].set_title("EPI")
		axs[2].set_title("ENDO")
		labels = ["EPI", "ENDO"]
		markers = ["+","x"]
		for n in 1:2
			axs[n].plot(samplePoint[1], samplePoint[2], marker=markers[n], color="k", markeredgewidth=1, markersize=3);
		end
		for n in 1:length(V)
			axs[3].plot(t, V[n], "-", label=labels[n], alpha=0.7)
		end
		for n in 1:length(V)
			axs[3].plot([t[sampleInt]], [V[n][sampleInt]], color="k", markeredgewidth=1, markersize=3, marker=markers[n])
		end
		labs = ["(a)", "(b)", "(c)"]
		for n in 1:3
			axs[n].text(-0.1, -0.1, labs[n], horizontalalignment="center", verticalalignment="center", transform=axs[n].transAxes)
		end
		xl = t[sampleInt].+[-550.0,+550.0];
		xl[1] = max(xl[1],t[1]); xl[2] = min(xl[2],t[end]);
		axs[3].set_xlim(xl)
		axs[3].set_ylim([-0.05, 1.2])
		axs[3].set_yticks([0.0,1.0])
		axs[3].legend(loc="upper center",edgecolor="none",ncol=2)
		axs[3].set_xlabel("\$ t \$ [ms]")
		plt.savefig("./$(ind)_datafig.pdf", bbox_inches="tight")
		plt.close()
	end
	
	# make a movie of the dynamics
	if movie
		recordingMovie(states, spaceMasks, timeMasks, mutualSpaceMask, mutualTimeMask; figname="$(ind)")
	end

	# output observations
	if outputObs
		# read in a single observation file as a template
		obs = readObsFile("./obs/template.dat")
		
		# make the uncertainties uniform from surfaces
		for n in 1:length(obs)
			obs[n][6] = 0.05
		end
				
		# extend obs by copying the elements with obs[n][2] == 1.0, with the rhs set to 4.0:3.0:50.0
		#=
		tmp = copy(obs[1])
		for z in 4.0:3.0:50.0, y in 1.0:3.0:200.0, x in 1.0:3.0:200.0
			push!(obs, Float32[tmp[1], z, y, x, 0.0, 0.0])
		end
		=#
		# extend obs by copying the elements with obs[n][2] == 1.0, with the rhs set to 24.0
		tmp = copy(obs[1]);
		obs = []
		for z in [1.0,50.0], y in 1.0:3.0:200.0, x in 1.0:3.0:200.0
			if (z != 24.0) err = 0.05 else err = 0.50 end
			push!(obs, Float32[tmp[1], z, y, x, 0.0, err])	# errors in obs array are passed through
		end
		
		# truncate the obs vector by those outside the respective masks
		#	Note: pass [mutualSpaceMask, mutualSpaceMask] as second argument to make observation domains the same
		checkObsAgainstMask!(obs, spaceMasks, corners)
		
		# generate the splines for the spatiotemporal states
		stateSplines = stateInterpolants(states, mutualTimeMask, mutualSpaceMask)
		
		# interpolate into stateSplines and mutate the observation template and save to new observation files 
		newObsDir = "./obs/$(ind)_lerp_dense/"; mkpath(newObsDir);
		sampleTimes = 2.0.*(1:length(mutualTimeMask))[mutualTimeMask] # 2.0 ms per discrete time index
		generateObservationFileSequence!(obs, stateSplines, newObsDir, "", sampleTimes, corners; err = 0.05, lerp = true)
	end
	
	return nothing
end

main(;plotting=false, outputObs=true, l=3.0)
