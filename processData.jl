#push!(LOAD_PATH,pwd())
include("loadAlessioData.jl")
using .loadAlessioData, DelimitedFiles, FortranFiles, Printf, ProgressBars, PyPlot, Peaks, Dierckx, Interpolations

function main(; plotting=true, movie=true)
	
	# where is the data
	basepath = "./2011-05-28_Rec78-103_Pace_Apex/"
	inds = getPossibleIndices(basepath)
	ind = rand(inds) # This is just useful for debugging; normally you would choose an index
	
	# load in some data 
	states = loadFile(basepath, ind)
	
	# form the space and time masks
	spaceMasks = [spatialMask(state) for state in states]
	timeMasks  = [temporalMask(state) for state in states]
	mutualSpaceMask = mutualMask(spaceMasks)
	mutualTimeMask = mutualMask(timeMasks)
	
	# print some diagnostic notices
	println("Loading data for $(ind)...")
	println("states stores EPI  data (recording, masks, etc) in first  index")
	println("states stores ENDO data (recording, masks, etc) in second index")
	
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
		plotState!(fig, axs, sampleTime, states, spaceMasks, timeMasks, [size(state) for state in states], mutualMask(spaceMasks), mutualMask(timeMasks); drawColorbar=true)
		plt.savefig("./state/$(ind)_state.svg", bbox_inches="tight")
		plt.close()
	end
	
	# make a movie of the dynamics
	if movie
		recordingMovie(states, spaceMasks, timeMasks; figname="$(ind)")
	end

	# read in observation files, change the obs / err values to interpolated values from input data, write to new obs file
	# need a loop here which takes 
	obs = readObsFile("./obs/0140.dat")
	stateSplines = stateInterpolants(states, mutualTimeMask, mutualSpaceMask)
	generateObservationFileSequence!(obs, stateSplines, "./obs/test/", "$(ind)", (1:length(mutualTimeMask)[mutualTimeMask]))
end

main()
