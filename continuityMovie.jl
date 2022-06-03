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

function continuityMovie(workDir, ens="gues", x=25, y=100, z=100; dx=0.015, dt=2.0)
	
	u = rand(Float64,nx,ny,nz)
	v = rand(Float64,nx,ny,nz)
	w = rand(Float64,nx,ny,nz)
	S = [view(u,:,:,:), view(v,:,:,:), view(w,:,:,:)]

	meanDir = "/data/$(workDir)/$(ens)/mean";
	frameDir = "./continuity_frames"; 
	mkpath(frameDir);
	
	nts = 1:length(readdir("$(meanDir)/")) 	# get how many gues time steps there are
	#=
	for n in 1:length(nts)				# check them by reading each corresponding anal dir
		try 
			readdir(@sprintf("%s/../sprd/%04d/",meanDir,n))
		catch
			nts = nts[1:(n-1)]
		end
	end
	=#	
	times = dt.*(nts.-1)
	
	cmaps = ["Oranges", "Purples", "Greens"]
	varnames = ["u", "v", "w"]

	fig,axs = plt.subplots(3,3, figsize=(SW,SW), sharex="col", sharey="row", constrained_layout=true)			
	for (n,t) in ProgressBar(enumerate(times))
		for ax in axs
			ax.cla()
			ax.axes.set_aspect("equal")
		end
		try
			fnames = [@sprintf("%s/%04d/restart3d.%03d",meanDir,n,p) for p in 0:3]
			readRestarts!(u, v, w, fnames)
			for o in 1:length(S)
				axs[o,3].pcolormesh(transpose(S[o][x,:,:]), snap=true, shading="auto", rasterized=true, vmin=0.0, vmax=1.0, cmap=cmaps[o])
				axs[o,2].pcolormesh(transpose(S[o][:,y,:]), snap=true, shading="auto", rasterized=true, vmin=0.0, vmax=1.0, cmap=cmaps[o])
				if n==1
					ima = axs[o,1].pcolormesh(transpose(S[o][:,:,z]), snap=true, shading="auto", rasterized=true, vmin=0.0, vmax=1.0, cmap=cmaps[o])
					clb = plt.colorbar(ima, ax=axs[o,:])
					clb.ax.set_title("\$ $(varnames[o])\$")
				else
					axs[o,1].pcolormesh(transpose(S[o][:,:,z]), snap=true, shading="auto", rasterized=true, vmin=0.0, vmax=1.0, cmap=cmaps[o])
				end
			end
			for ax in axs
				ax.set_xticks([])
				ax.set_yticks([])
			end
			plt.suptitle("\$t=$(t)\$ [ms]")
			axs[1,3].set_title("\$$(x)/$(nx)\$")
			axs[1,2].set_title("\$$(y)/$(ny)\$")
			axs[1,1].set_title("\$$(z)/$(nz)\$")
			plt.savefig(@sprintf("%s/%04d.png", frameDir, n), dpi=300, bbox_inches="tight")
		catch
			print("\n\t Stopping frame generation at iteration $(n).\n")
			break
		end
	end
	try
		rm("$(workDir)_continuity.mp4");
	catch
	end
	run(`ffmpeg -r 10 -f image2 -start_number 1 -i $(frameDir)/%04d.png -c:v h264_nvenc -pix_fmt yuv420p -loglevel quiet -stats $(workDir)_continuity.mp4`);
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
	#print("This is >1.5-1.3× faster than the Python implementation!\n")	
	continuityMovie(workDir)    
	print("\n")
end

main()
