include("loadAlessioData.jl")
using .loadAlessioData, DelimitedFiles, FortranFiles, Printf, ProgressBars, PyPlot, Peaks, Dierckx, Interpolations

function main(; plotting=false, movie=false, outputObs=true)
	
	# where is the data
	basepath = "./2011-05-28_Rec78-103_Pace_Apex/"
	inds = getPossibleIndices(basepath)
	ind = rand(inds) # This is just useful for debugging; normally you would choose an index
	
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
	corners = loadAlessioData.square(theta=0.0, l=3.0, center=center)
	
	# make some plots
	if plotting
	# sampling position is randomly sampled and checked against mutual space mask
		samplePoint = [rand(1:128), rand(1:128)]
		while mutualSpaceMask[samplePoint[1],samplePoint[2]] == false
			samplePoint = [rand(1:128), rand(1:128)]
		end
		t, V = sampleData(states, samplePoint, mutualTimeMask, mutualSpaceMask)
		fig, axs = plt.subplots(1,3,figsize=(12,3),constrained_layout=true)
		plotTrace!(fig, axs, t, V; threshold=0.25)
		plt.savefig("./trace/$(ind)_trace.svg", bbox_inches="tight")
		plt.close()
	
		# plotting a single state at a prescribed time
		sampleTime = rand((1:length(mutualTimeMask))[mutualTimeMask])
		fig, axs = plt.subplots(1,2,figsize=(8,3),constrained_layout=true)
		plotState!(fig, axs, sampleTime, states, spaceMasks, timeMasks, [size(state) for state in states], mutualSpaceMask, mutualTimeMask; drawColorbar=true)
		plotSamplingSquare!(fig, axs, corners)
		plt.savefig("./state/$(ind)_state.svg", bbox_inches="tight")
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
		
		# generate the splines for the spatiotemporal states
		stateSplines = stateInterpolants(states, mutualTimeMask, mutualSpaceMask)
		
		# interpolate into stateSplines and mutate the observation template and save to new observation files 
		newObsDir = "./obs/$(ind)/"; mkpath(newObsDir);
		sampleTimes = 2.0.*(1:length(mutualTimeMask))[mutualTimeMask]
		generateObservationFileSequence!(obs, stateSplines, newObsDir, "", sampleTimes; center = center)
	end
	
	return nothing
end

main(;plotting=true, outputObs=false)
