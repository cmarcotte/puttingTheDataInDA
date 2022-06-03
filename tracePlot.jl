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

function tracePlot(workDir, obsDir; zs=[1],ys=[160],xs=[100], dx=0.015, dt=2.0)

	gu = zeros(Float64,nx,ny,nz)
	gv = zeros(Float64,nx,ny,nz)
	gw = zeros(Float64,nx,ny,nz)
	au = zeros(Float64,nx,ny,nz)
	av = zeros(Float64,nx,ny,nz)
	aw = zeros(Float64,nx,ny,nz)
	
	guesDir = "/data/$(workDir)/gues/mean";
	analDir = "/data/$(workDir)/anal/mean";	
	nts = 1:length(readdir("$(guesDir)/")) 	# get how many gues time steps there are
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
	@assert length(zs) == length(ys) == length(xs)
	trcs = zeros(Float64, length(times), 7, length(zs))
	
	for (n,t) in ProgressBar(enumerate(times))
		fnames = [@sprintf("%s/%04d/restart3d.%03d",guesDir,n,p) for p in 0:3]
		readRestarts!(gu, gv, gw, fnames)
		fnames = [@sprintf("%s/%04d/restart3d.%03d",analDir,n,p) for p in 0:3]
		readRestarts!(au, av, aw, fnames)
		obs = reduce(hcat,readObsFile(@sprintf("%s/%04d.dat",obsDir,n)))
		
		for o in eachindex(zs)
			
			trcs[n,1,o] = gu[zs[o],ys[o],xs[o]]
			trcs[n,2,o] = gv[zs[o],ys[o],xs[o]]
			trcs[n,3,o] = gw[zs[o],ys[o],xs[o]]
			
			trcs[n,4,o] = au[zs[o],ys[o],xs[o]]
			trcs[n,5,o] = av[zs[o],ys[o],xs[o]]
			trcs[n,6,o] = aw[zs[o],ys[o],xs[o]]
			
			iz = findall(obs[2,:].==zs[o] .&& obs[3,:].==ys[o] .&& obs[4,:].==xs[o])[]
			trcs[n,7,o] = obs[5,iz]
		end
	end
	
	fig, axs = plt.subplots(3, 1, figsize=(SW, SW), sharex=true, sharey=true, constrained_layout=true)
	labs = ["(a)", "(b)", "(c)"]
	vrns = ["\$u\$", "\$v\$", "\$w\$"]
	for o in 1:length(zs)
		for m in 1:3
			axs[m,o].plot(times, trcs[:,m,o],   "-C0", label="Background")
			axs[m,o].plot(times, trcs[:,3+m,o], "-C1", label="Analysis")
			if m==1
				axs[m,o].plot(times, trcs[:,end,o], ".-C2", linewidth=0.5, markersize=4, label="Obs. $([zs[o],ys[o],xs[o]])")
			end
			axs[m,o].set_ylabel("$(vrns[m])")
			axs[m,o].text(-0.2, +1.1, labs[m], horizontalalignment="center", verticalalignment="center", transform=axs[m,o].transAxes)
		end
		axs[1,o].legend(loc="lower left", bbox_to_anchor= (0.0, 1.01), ncol=1, borderaxespad=0, frameon=false)
		#axs[1,o].set_title("\$[x,y,z]=$([zs[o],ys[o],xs[o]])\$")
		axs[end,o].set_xlabel("\$t\$ [ms]")
		axs[end,o].set_xlim([0,2000])
	end
	plt.savefig("$(workDir)_traces.pdf", bbox_inches="tight", dpi=600)
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
	
	#print("This is >3.33-1.05× faster than the Python implementation!\n")	
	tracePlot(workDir, obsDir)    
	print("\n")
end

main()
