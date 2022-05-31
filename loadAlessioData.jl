module loadAlessioData
#=
	This module implements convenience functions for
	- loading of Alessio's data
	- reading and writing Fortran observation data files
	- plotting and analysis of the data files 
	- interpolating Alessio's data onto appropriate grids
=#

# For:	     data input/output             Printing        Plotting        Data Analysis
#      /--------------------------\  /------------------\  /----\  /----------------------------\   
using DelimitedFiles, FortranFiles, Printf, ProgressBars, PyPlot, Peaks, Dierckx, Interpolations

export getPossibleIndices, loadFile, spatialMask, temporalMask, mutualMask, sampleData, plotTrace!, plotState!, recordingMovie, readObsFile, stateInterpolants, generateObservationFileSequence!, spatialMaskCentroid, plotSamplingSquare!

#
# Convenience functions to identify the available data files and load the data consistently
#
function getPossibleIndices(basepath)
	# lists possible data file indices
	files = readdir(basepath)
	indices = []
	for file in files
		if cmp(file[1:6], "Erased") == 0
			push!(indices, parse(Int,file[29:31]))
		end
	end
	sort!(indices)
	unique!(indices)
	
	return indices
end

function loadData(filename) # loadData function code due to Max Comstock
	# Pre-allocate data array
	total_elements = Int(stat(filename).size / sizeof(UInt16))
	obs_data = Array{UInt16, 1}(undef, total_elements)
	
	# Read the data and convert from big-endian to system endianness
	read!(filename, obs_data)
	obs_data = ntoh.(obs_data)
	
	# Shape into frames of 128×128 readings
	obs_data = reshape(obs_data, :, 128, 128)

	# rescale data to [0,1] and convert to Float64 (overkill; Float32 would work, Float16 is too small) 
	obs_data = obs_data ./ (maximum(obs_data[:])-minimum(obs_data[:]))
	
	return obs_data
end

function transformData!(obs)
	# transforms data in-place for consistent orientation; array is now [EPI,ENDO][t,x,y]
	obs[1] .= obs[1][:,:,end:-1:1]
	obs[2] .= obs[2][:,end:-1:1,:]
	return nothing
end

function loadFile(basepath, index::Int)
	# loads data from a pair of files
	basefile = "Erased_2011-05-18_Exp000_Rec"*@sprintf("%03i",index)
	# EPI first, then ENDO
	filenames= [basepath*basefile*"_Cam0-Yellow", basepath*basefile*"_Cam1-Blue"]
	dats = []
	for fname in filenames
		dat = loadData(fname)
		push!(dats, dat)
	end
	
	# re-orient data arrays for consistency
	transformData!(dats)
	
	return dats
end

#
# Basic analysis functions to generate masks for the experimental data arrays
#

function spatialMask(obs_data)
	# generates BitArray of size [Nx,Ny] which corresponds to valid spatial indices
	N = size(obs_data)
	mask = zeros(Bool, N[2], N[3])
	for n in 1:N[2], m in 1:N[3]
		for t in 1:N[1]
			if mask[n,m] == false && obs_data[t,n,m] > 0
				mask[n,m] = 1
			end
		end
	end
	
	return mask		
end

function temporalMask(obs_data)
	# generates BitArray of size [Nt] which corresponds to valid temporal indices
	N = size(obs_data)
	mask = zeros(Bool, N[1])
	for t in 1:N[1]
		for n in 1:N[2], m in 1:N[3]
			if mask[t] == false && obs_data[t,n,m] > 0
				mask[t] = 1
			end
		end
	end
	
	return mask		
end

function mutualMask(mask1, mask2)
	# pairwise function to reduce the data ranges for the EPI and ENDO recordings to a consistent subset 
	if length(size(mask1)) == 1 && length(size(mask2)) == 1 && length(mask1) != length(mask2)
		println("Temporal masks are of different lengths; truncating to shorter")
		tmp = min(length(mask1), length(mask2))
		mask1 = mask1[1:tmp]
		mask2 = mask2[1:tmp]
	end
	@assert size(mask1) == size(mask2)
	
	return mask1 .* mask2
end

function mutualMask(masks)
	# array input of masks
	n = length(masks)
	mutmask = masks[1]
	for m = 2:n
		mutmask = mutualMask(mutmask, masks[m])
	end
	
	return mutmask 
end

function spatialMaskCentroid(mask; dx=0.06)

	N = size(mask)
	ctr = 0.0;
	for x in 1:N[1], y in 1:N[2]
		ctr += (x + y*im)*mask[x,y];
	end
	return Int.(round.([imag(ctr), real(ctr)]./sum(mask))).*dx # note the swap in real/imag
end

#
# Sampling functions
#

function sampleData(obs, sampleIndex::Array{T,1}, mutualTimeMask, mutualSpaceMask; dt=2.0) where T <: Integer
	# generates (t,V(t)) arrays for EPI and ENDO recordings in obs
	# first check that the spatial indices are valid
	@assert mutualSpaceMask[sampleIndex[1],sampleIndex[2]] == true
	# generate times from the mutualTimeMask
	t = (1:length(mutualTimeMask))[mutualTimeMask]
	# select obs[t, sampleIndex[1], sampleIndex[2]] for EPI and ENDO
	V = [obs[n][t, sampleIndex[1], sampleIndex[2]] for n=1:2]
	
	return (t.*dt, V)
end

function sampleData(obs, samplePoint::Array{T,1}, mutualTimeMask, mutualSpaceMask; h=0.02) where T <: AbstractFloat
	sampleIndex = Int.(round.(realSamplePoint ./ h))
	return sampleData(obs, sampleIndex, mutualTimeMask, mutualSpaceMask)
end

#
# Analysis functions for the data
#

function analyzeVoltage(t, V; threshold=0.25)
	#=
	This function uses Dierckx.jl for interpolations and root-finding. 
	I would like to use Interpolations.jl, but there's no in-built root-
	finding, which would mean using Roots.jl, so might as well just use Dierkx.jl.
	=#
	# Interpolate the voltage (at a specific spatial point) over time 
	Vt = Spline1D(t[:], V.-threshold, k=3);
	R = roots(Vt; maxn=Int(5e3));		# time points R: V(R) == threshold
	D = derivative(Vt, R);		# V'(R)

	# storage for APD, DI, APA for this V(t)
	APD = Float64[]
	DI  = Float64[]
	APA = Float64[]
	
	# for each root, check whether V'(R)≷ 0 to detect beginning or end of AP
	for n in 1:length(R)-1
	       if D[n] > 0 && D[n+1] < 0
		       push!(APD, R[n+1]-R[n])
		       push!(APA, threshold+maximum(Vt(R[n]:R[n+1])))
	       elseif D[n] < 0 && D[n+1] > 0
		       push!(DI, R[n+1]-R[n])
	       end
	end
	return (APD, DI, APA)
end

function peaks(V; minimumProminence=0.9)

	pks, vals = findmaxima(V)
	pks, proms = peakproms(pks, V)
	
	inds = (1:length(pks))[proms .>= minimumProminence]
	
	return (inds, pks, proms)
end

# this is convenient, but expensive, and the interpolant throws when (t,x,y,z) are outside the available ranges, but doesn't tell you what those ranges are.
function generateCompleteInterpolant(states, mutualSpaceMask, mutualTimeMask; dt=2.0, dx=0.06, dz=1.0) # ms and cm, respectively; dz is thickness of tissue
	# generates a single data interpolant across time and space)
	it = (1:length(mutualTimeMask))[mutualTimeMask]
	u = cat(states[1][it,:,:], states[2][it,:,:]; dims=4) # huge temporary array; probably a way to avoid?
	itp = LinearInterpolation((dt.*(0:length(it)-1), (0:127).*dx, (0:127).*dx, [0,dz]), u) # repackages that same huge array
	# now itp(t::T, x::T, y::T, z::T) where T <: Real gives u(t,x,y,z) for the dataset
	# we we'd like to filter out calls that are outside mutualSpaceMask (the mutualTimeMask already throws in the interpolant).
	# so we can first construct a function which filters the space indices as well.
	#function interpolant(t,x,y,z)
	#	if mutualSpaceMask[Int(round(y/dx)),Int(round(x/dx))]
	#		return itp(t,x,y,z)
	#	else
	#		error("(x,y) coordinate outside mutualSpaceMask")
	#	end	
	#end
	#return interpolant
	return itp
end

function stateInterpolants(states, mutualTimeMask, mutualSpaceMask; dt=2.0, dx=0.06, dy=0.06)
	
	tlims = (1:length(mutualTimeMask))[mutualTimeMask][[1;end]]
	t0 = dt.*(tlims[1]:tlims[2])
	x0 = dx.*((1:size(mutualSpaceMask,1)).-1)
	y0 = dy.*((1:size(mutualSpaceMask,2)).-1)
	
	stateSplines = [ CubicSplineInterpolation((t0, x0, y0), state[(1:length(mutualTimeMask))[mutualTimeMask],:,:]) for state in states ]
	
	return stateSplines
end

function square(;theta=0.0, l=3.0, center=[3.5;3.5])
	# generates the bounding corners (with first and last repeating) of the square for sampling the states
	if theta ≠ 0.0
		@warn "θ ≠ 0 is currently unimplemented, sorry!"
	end
	if abs(theta) >= pi/2
		@warn "θ is mapped to (-π/2,+π/2)"
		theta = mod(theta, pi/2)
	end
	corners = center .+ l .* ([cos(theta) sin(theta); -sin(theta) cos(theta)] * [-1 1 1 -1; -1 -1 1 1]./2);
	corners = hcat(corners[:,:], corners[:,1])

	return corners
end

function coordinateVectors(corners)
	# corners should be generated from the square function
	# corners = square(theta=theta, l=l, center=center)
	x⃗ = corners[:,2].-corners[:,1]
	y⃗ = corners[:,4].-corners[:,1]
	# s.t.  ∀(a,b) ∈ [-0.5, +0.5]: center .+ a.*x⃗ .+ b.*y⃗ ∈ square
	return (x⃗, y⃗)
end

#
# Data plotting functions
#

function plotTrace!(fig, axs, t, V; threshold=0.5) where T <: Integer

	labels = ["EPI", "ENDO"]
	markers = ["+","x"]
	
	for ax in axs
		ax.cla()
	end
	
	for n in 1:length(V) # [EPI, ENDO]
		axs[1].plot(t, V[n], "-", label=labels[n], alpha=0.7)
		inds, pks, proms = peaks(V[n]; minimumProminence=0.1)
		axs[1].plot(t[pks[inds]], V[n][pks[inds]], markers[n]*"k")
		axs[1].plot(t, threshold.*ones(size(t)), "--k", label="", alpha=0.7)
		
		BCL = diff(t[pks[inds]]); 
		binwidth=2.0; bins=((round(minimum(BCL)/binwidth)-1)*binwidth):binwidth:((round(maximum(BCL)/binwidth)+1)*binwidth);
		axs[2].hist(BCL, bins=bins, density=true, log=true, label=labels[n], color="C$(n-1)", alpha=0.7)
		
		APD, DI, APA = analyzeVoltage(t, V[n]; threshold=threshold)
		axs[3].errorbar(APD, APA, xerr=2.0, yerr=0.05, fmt=".", label=labels[n], alpha=0.7)
		
		cBCL = min(length(APD),length(DI)); 
		BCL = APD[1:cBCL] .+ DI[1:cBCL]
		axs[2].hist(BCL, bins=bins, density=true, histtype="step", color="C$(n-1)", log=true, label="APD+DI,"*labels[n], alpha=0.7)
	end
	axs[1].set_xlabel("\$ t \$ [ms]")
	axs[1].set_ylabel("\$ V(t) \$ [a.u.]")
	axs[1].legend(loc=0,edgecolor="none")
	axs[2].set_xlabel("BCL [ms]")
	axs[2].set_ylabel("\$P(\$BCL\$)\$")
	axs[2].legend(loc=0,edgecolor="none", ncol=2)
	axs[3].set_xlabel("APD [ms]")
	axs[3].set_ylabel("APA [a.u.]")
	axs[3].legend(loc=0,edgecolor="none")
	
	return nothing
end

function plotState!(fig, axs, t, obs, spaceMasks, timeMasks, N, spacemask, timemask; drawColorbar=false)
	
	if drawColorbar  # only draw the colorbar on the first frame
		imm = axs[1].imshow(obs[1][t,:,:], aspect="equal", vmin=0.0, vmax=1.0, cmap="Oranges", interpolation="none", origin="lower")	
		clb = fig.colorbar(imm, ax=axs, aspect=50)
	end
	labels = ["EPI", "ENDO"]
	
	for n=1:2
		axs[n].cla()
		axs[n].imshow(obs[n][t,:,:], alpha=Float64.(spaceMasks[n]), aspect="equal", vmin=0.0, vmax=1.0, cmap="Oranges", interpolation="none", origin="lower")
		axs[n].contour(Float64.(spacemask), "k", levels=[0.5], linewidths=0.5)
		axs[n].axis("off")
		axs[n].set_xticks([])
		axs[n].set_yticks([])
		axs[n].set_title("$(labels[n]), \$ t=\$"*"$(t*2) ms")
	end
	axs[1].contour(Float64.(spaceMasks[2]), "k", levels=[0.5], linestyles="dashed", linewidths=0.5)
	axs[2].contour(Float64.(spaceMasks[1]), "k", levels=[0.5], linestyles="dashed", linewidths=0.5)

	return nothing	
end

function plotSamplingSquare!(fig, axs, corners)

	axs[1].plot(corners[1,:]./0.06, corners[2,:]./0.06, "--k", linewidth=1)
	axs[2].plot(corners[1,:]./0.06, corners[2,:]./0.06, "--k", linewidth=1)

	return nothing
end

function recordingMovie(obs, spaceMasks, timeMasks, spacemask, timemask; figname="recording")
	
	N = size(obs[1])
	
	rm("./frames/", force=true, recursive=true)
	mkdir("./frames/")
	
	println("Making movie frames:")
	fig,axs = PyPlot.subplots(1,2,figsize=(8,3),constrained_layout=true)
	
	for t in ProgressBar((1:length(timemask))[timemask])
		
		plotState!(fig, axs, t, obs, spaceMasks, timeMasks, N, spacemask, timemask; drawColorbar=(t==(1:length(timemask))[timemask][1]))
		
		imname = @sprintf("%04i",t)
		PyPlot.savefig("./frames/$(figname)_$(imname).png")
	end
	
	PyPlot.close("all")
	
	println("Running frames through FFMPEG:")
	t = (1:length(timemask))[timemask][1]
	ffmpegCommand = `ffmpeg -r 60 -f image2 -start_number $(t) -i ./frames/$(figname)_%04d.png -c:v hevc -crf 25 -pix_fmt yuv420p ./dynamics/$(figname)_dynamics.mp4`;
	println(ffmpegCommand)
	run(ffmpegCommand);
	
	rm("./frames/", force=true, recursive=true)
	
	return nothing
end

#
# Observation file functions
#

function readObsFile(obsFile::String; verbose=false)
	#=
	obsFile : string for observation filename
	=#
	f = FortranFile(obsFile)
	#obs = Array{Float32,1}(undef, 6) # not needed, but useful for reference
	Obs = [] # accumulator
	while !eof(f)
		obs = read(f, (Float32,6));
		if verbose; println("Read observation $(obs)"); end
		push!(Obs, obs)
	end
	if eof(f)
		close(f)
	end
	return Obs
end

function fillObservations!(obs, stateSplines, t; theta=0.0, l=3.0, center=[3.5;3.5], min_x = 0, max_x = 200, min_y = 0, max_y = 200, err = 0.025)
	# states is a Array{Array[Float64,3},1}
	#	states = [EPI, ENDO]
	# Obs is a Array{Array{Float32,1},1} where each element of the top array is of form:
	#	[id, x, y, z, val, err] = [2819.0   1.00  89.00    100.00    2.15E-01    5.00E-02]
	# this updates, in-place, the val and err according to the interpolated state at (x,y,z)
	
	# the data has shape nt x 2 x 128 x 128
	# while the observations are 1 x 2 x 200 x 200
	# so we need to interpolate the data into two continuous (x,y) functions
	# and evaluate the functions at
	# the interpolated relative observation positions {0,1} x [0,1], x [0,1]
	# where the data has (0,0), (0,1), (1,0), and (1,1) at a prescribed and inscribed square
	
	corners = square(theta=theta, l=l, center=center)

	for n in 1:length(obs)
		
		# get relative position (within square)
		#x = corners[1,1] + l*cos(theta)*(obs[n][3]-min_x)/(max_x-min_x)
		#y = corners[2,1] + l*sin(theta)*(obs[n][4]-min_y)/(max_y-min_y)
		x = corners[1,1] + (corners[1,2]-corners[1,1])*(obs[n][3]-min_x)/(max_x-min_x)
		y = corners[2,1] + (corners[2,3]-corners[2,1])*(obs[n][4]-min_y)/(max_y-min_y)
		
		if obs[n][2] == 1.0		# if on bottom => EPI sample
			obs[n][5] = stateSplines[1](t,y,x)
		elseif obs[n][2] == 50.0	# if on top => ENDO sample
			obs[n][5] = stateSplines[2](t,y,x)
		end
		
		obs[n][6] = err 		# What is the error in the observation here
	end
	return nothing				# obs is modified in-place
end

function writeObsFile(obsFile::String, Obs)
	#=
	obsFile : string for observation filename
	Obs: observation array
	=#
	f = FortranFile(obsFile, "w")
	#obs = Array{Float32,1}(undef, 6) # not needed, but useful for reference
	for n in 1:length(Obs)
		write(f, Obs[n])
		@debug "Wrote observation $(Obs[n])"
	end
	close(f)
	@debug "Wrote observation file $(obsFile)"
	return nothing
end

function generateObservationFileSequence!(obs, stateSplines, obsDir, sequenceName, sampleTimes; theta=0.0, l=3.0, center=[3.5;3.5], min_x = 0, max_x = 200, min_y = 0, max_y = 200, err = 0.025)
	#=
	We loop over the sampleTimes, modifying the template observation array obs, and writing the array to a new obs file	
	=#
	@show center
	println("Writing new observation files to $(obsDir):")
	for (n,t) in ProgressBar(enumerate(sampleTimes))
		fillObservations!(obs, stateSplines, t; theta=theta, l=l, center=center, min_x = min_x, max_x = max_x, min_y = min_y, max_y = max_y, err = 0.025)
		obsFile = @sprintf("%s/%s%04d.dat", obsDir, sequenceName, n)
		writeObsFile(obsFile, obs)
	end
	return nothing
end

end
