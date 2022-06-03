include("readData.jl")
using .readData, FortranFiles, PyPlot, Printf, ArgParse, ProgressBars, SparseArrays, Statistics

plt.style.use("seaborn-paper")

PyPlot.rc("font", family="serif")
PyPlot.rc("text", usetex=true)
PyPlot.matplotlib.rcParams["axes.titlesize"] = 10
PyPlot.matplotlib.rcParams["axes.labelsize"] = 10
PyPlot.matplotlib.rcParams["xtick.labelsize"] = 9
PyPlot.matplotlib.rcParams["ytick.labelsize"] = 9

const SW = 3.40457
const DW = 7.05826

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

function obs2arr(obs)
	return reduce(hcat,obs)
end

function errorPlot(workDir, obsDir; zs=[1,50], ens="gues", dx=0.015, dt=2.0)

	mu = zeros(Float64,nx,ny,nz)
	mv = zeros(Float64,nx,ny,nz)
	mw = zeros(Float64,nx,ny,nz)
	su = zeros(Float64,nx,ny,nz)
	sv = zeros(Float64,nx,ny,nz)
	sw = zeros(Float64,nx,ny,nz)
	
	meanDir = "/data/$(workDir)/$(ens)/mean";
	sprdDir = "/data/$(workDir)/$(ens)/sprd";	
	nts = 1:length(readdir("$(meanDir)/")) 	# get how many gues time steps there are
	#=
	for n in 1:length(nts)				# check them by reading each corresponding anal dir
		try 
			readdir(@sprintf("%s/%04d/",sprdDir,n))
		catch
			nts = nts[1:(n-1)]
		end
	end
	=#	
	times = dt.*(nts.-1)
	
	# accumulation arrays
	errrms = zeros(Float64, length(times), length(zs))
	errstd = zeros(Float64, length(times), length(zs))
	sprrms = zeros(Float64, length(times), length(zs))
	sprstd = zeros(Float64, length(times), length(zs))
	
	for (n,t) in ProgressBar(enumerate(times))
		fnames = [@sprintf("%s/%04d/restart3d.%03d",meanDir,n,p) for p in 0:3]
		readRestarts!(mu, mv, mw, fnames)
		fnames = [@sprintf("%s/%04d/restart3d.%03d",sprdDir,n,p) for p in 0:3]
		readRestarts!(su, sv, sw, fnames)
		obs = reduce(hcat,readObsFile(@sprintf("%s/%04d.dat",obsDir,n)))
		
		for (o,z) in enumerate(zs)
			iz = findall(obs[2,:].==zs[o]) #obs[2,iz].==zs[o], so obs[5,iz] are the u-values
			X = Int.(obs[2,iz]); Y = Int.(obs[3,iz]); Z = Int.(obs[4,iz]); U = obs[5,iz];
			tmp = zeros(Float64, length(iz), 2)
			for l in eachindex(iz)
				tmp[l,1] = mu[X[l],Y[l],Z[l]]-U[l]
				tmp[l,2] = su[X[l],Y[l],Z[l]]
			end
			errrms[n,o] = sqrt.(sum(abs2, tmp[:,1])/length(iz))
			errstd[n,o] = std(tmp[:,1])
			sprrms[n,o] = sqrt.(sum(abs2, tmp[:,2])/length(iz))
			sprstd[n,o] = std(tmp[:,2])
		end
	end

	fig, axs = plt.subplots(2, 1, figsize=(SW, SW), sharex=true, sharey=true, constrained_layout=true)
	labels = ["EPI", "ENDO"]
	labs = ["(a)", "(b)"]
	for (o,z) in enumerate(zs)
		axs[o].fill_between(times,  errrms[:,o] .- errstd[:,o],  errrms[:,o] .+ errstd[:,o], color="C0", linewidth=0, alpha=0.3, label="")
		axs[o].fill_between(times,  sprrms[:,o] .- sprstd[:,o],  sprrms[:,o] .+ sprstd[:,o], color="C1", linewidth=0, alpha=0.3, label="")
		axs[o].plot(times, errrms[:,o],  "-C0", label="RMS Error")
		axs[o].plot(times, sprrms[:,o], "-C1", linewidth=1, label="RMS Spread")
		axs[o].set_ylabel("$(labels[o])")
		axs[o].text(-0.2, +1.1, labs[o], horizontalalignment="center", verticalalignment="center", transform=axs[o].transAxes)
	end
	axs[begin].legend(loc="lower left", bbox_to_anchor= (0.0, 1.01), ncol=2, borderaxespad=0, frameon=false)
	axs[end].set_xlabel("\$t\$ [ms]")
	axs[end].set_xlim([0,2000])
	axs[begin].set_ylim([0.0,1.0])
	plt.savefig("$(workDir)_SurfaceErrors.pdf", bbox_inches="tight", dpi=300)
	plt.close()
	
	return nothing
end

function parse_commandline()
	s = ArgParseSettings()

	@add_arg_table s begin
	"workDir"
		help = "workDir"
		required = true
	end

	return parse_args(s)
end

function main(; baseDir="/data")
	parsed_args = parse_commandline()
	workDir = parsed_args["workDir"]
	obsDir = "/data/$(workDir)/obs"
	
	#print("This is >3.86-1.17× faster than the Python implementation!\n")	
	errorPlot(workDir, obsDir)    
	print("\n")
end

main()
