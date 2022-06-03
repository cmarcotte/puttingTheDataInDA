include("readData.jl")
using Main.readData, PyPlot, Printf, ArgParse, ProgressBars, FortranFiles

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

function dynamicsMovie(workDir, obsDir, ens = "gues", zs = [1,24,50]; dx=0.015, dt=2.0)
		
	mu = rand(Float64,nx,ny,nz)
	mv = rand(Float64,nx,ny,nz)
	mw = rand(Float64,nx,ny,nz)
	su = rand(Float64,nx,ny,nz)
	sv = rand(Float64,nx,ny,nz)
	sw = rand(Float64,nx,ny,nz)
	
	meanDir = "/data/$(workDir)/$(ens)/mean";
	sprdDir = "/data/$(workDir)/$(ens)/sprd";
	frameDir = "./dynamics_frames"; 
	mkpath(frameDir);
	
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
	
	fig,axs = plt.subplots(3,3, figsize=(SW,SW), sharex=true, sharey=true, constrained_layout=true)			
	for (n,t) in ProgressBar(enumerate(times))
		for ax in axs
			ax.cla()
			ax.axes.set_aspect("equal")
		end
		try
			fnames = [@sprintf("%s/%04d/restart3d.%03d",meanDir,n,p) for p in 0:3]
			readRestarts!(mu, mv, mw, fnames)
			fnames = [@sprintf("%s/%04d/restart3d.%03d",sprdDir,n,p) for p in 0:3]
			readRestarts!(su, sv, sw, fnames)
			obs = reduce(hcat,readObsFile(@sprintf("%s/%04d.dat",obsDir,n)))
			
			for (o,z) in enumerate(zs)
				obsdepth = z==1 ? 0.0 : z*dx
				axs[end,o].set_xlabel(@sprintf("\$ %2.3f\$ [cm]", obsdepth))
				for m in 1:3
					axs[m,o].set_xticks([])
					axs[m,o].set_yticks([])
				end
				axs[o,1].pcolormesh(transpose(mu[z,:,:]), snap=true, shading="auto", rasterized=true, vmin=0.0, vmax=1.0, cmap="Oranges")
				axs[o,2].pcolormesh(transpose(su[z,:,:]), snap=true, shading="auto", rasterized=true, vmin=0.0, vmax=1.0, cmap="Oranges")
				iz = findall(z.==obs[2,:])
				if length(iz) > 0
					axs[o,3].tricontourf(obs[3,iz][:], obs[4,iz][:], obs[5,iz][:], levels=129, vmin=0.0, vmax=1.0, cmap="Oranges")
				end
			end
			if n==1
				ima = axs[length(zs),1].pcolormesh(transpose(mu[zs[end],:,:]), snap=true, shading="auto", rasterized=true, vmin=0.0, vmax=1.0, cmap="Oranges")
				clb = plt.colorbar(ima, ax=axs)
				clb.ax.set_title("\$u\$")
				axs[1,1].set_xlim([1,200])
				axs[1,1].set_ylim([1,200])
			end
			axs[1,1].set_title("\$ \$Mean")
			axs[1,2].set_title("\$ \$Spread")
			axs[1,3].set_title("\$ \$Obs")
			plt.suptitle("\$t=$(t)\$ [ms]")
			plt.savefig(@sprintf("%s/%04d.png", frameDir, n), dpi=300, bbox_inches="tight")
		catch
			print("\n\t Stopping frame generation at iteration $(n).\n")
			break
		end
	end
	try
		rm("$(workDir)_dynamics.mp4");
	catch
	end
	run(`ffmpeg -r 10 -f image2 -start_number 1 -i $(frameDir)/%04d.png -c:v h264_nvenc -pix_fmt yuv420p -loglevel quiet -stats $(workDir)_dynamics.mp4`);
	rm(frameDir, recursive=true);
	return nothing
end

function parse_commandline()
	s = ArgParseSettings()

	@add_arg_table s begin
	"workDir"
		help = "workDir"
		required = true
	#"obsDir"
	#	help = "obsDir"
	#	required = false
	end

	return parse_args(s)
end

function main(; baseDir="/data")
	parsed_args = parse_commandline()
	workDir = parsed_args["workDir"]
	#try
	#	obsDir = parsed_args["obsDir"]
	#catch
		obsDir = "/data/$(workDir)/obs"
	#end
	#print("This is >5.43-1.22× faster than the Python implementation!\n")	
	dynamicsMovie(workDir, obsDir)    
	print("\n")
end

main()
