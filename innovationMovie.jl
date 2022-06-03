include("readData.jl")
using Main.readData, PyPlot, Printf, ArgParse, ProgressBars

plt.style.use("seaborn-paper")

PyPlot.rc("font", family="serif")
PyPlot.rc("text", usetex=true)
PyPlot.matplotlib.rcParams["axes.titlesize"] = 10
PyPlot.matplotlib.rcParams["axes.labelsize"] = 10
PyPlot.matplotlib.rcParams["xtick.labelsize"] = 9
PyPlot.matplotlib.rcParams["ytick.labelsize"] = 9

const sw = 3.40457
const dw = 7.05826

function innovationMovie(workDir, zs = [1,10,20,30,40,50]; dx=0.015, dt=2.0)
		
	gu = rand(Float64,nx,ny,nz)
	gv = rand(Float64,nx,ny,nz)
	gw = rand(Float64,nx,ny,nz)
	au = rand(Float64,nx,ny,nz)
	av = rand(Float64,nx,ny,nz)
	aw = rand(Float64,nx,ny,nz)
	
	guesDir = "/data/$(workDir)/gues";
	analDir = "/data/$(workDir)/anal";
	frameDir = "./innovation_frames"; 
	mkpath(frameDir);
	
	nts = 1:length(readdir("$(guesDir)/mean/")) 	# get how many gues time steps there are
	#=
	for n in 1:length(nts)				# check them by reading each corresponding anal dir
		try 
			readdir(@sprintf("%s/mean/%04d/",analDir,n))
		catch
			nts = nts[1:(n-1)]
		end
	end
	=#	
	times = dt.*(nts.-1)
	
	fig,axs = plt.subplots(3,length(zs), figsize=(dw,sw), sharex=true, sharey=true, constrained_layout=true)			
	for (n,t) in ProgressBar(enumerate(times))
		for ax in axs
			ax.cla()
			ax.axes.set_aspect("equal")
		end
		try
			fnames = [@sprintf("%s/mean/%04d/restart3d.%03d",guesDir,n,p) for p in 0:3]
			readRestarts!(gu, gv, gw, fnames)
			fnames = [@sprintf("%s/mean/%04d/restart3d.%03d",analDir,n,p) for p in 0:3]
			readRestarts!(au, av, aw, fnames)
			
			for (o,z) in enumerate(zs)
				obsdepth = z==1 ? 0.0 : z*dx
				axs[end,o].set_xlabel(@sprintf("\$ %2.3f\$ [cm]", obsdepth))
				for m in 1:3
					axs[m,o].set_xticks([])
					axs[m,o].set_yticks([])
				end
				if n==1 && o==length(zs)
					ima = axs[1,o].pcolormesh(transpose(gu[z,:,:]), snap=true, shading="auto", rasterized=true, vmin=0.0, vmax=1.0, cmap="Oranges")
					clb = plt.colorbar(ima, ax=axs[1,:])
					clb.ax.set_ylabel("\$u^b\$")
					ima = axs[2,o].pcolormesh(transpose(au[z,:,:]-gu[z,:,:]), snap=true, shading="auto", rasterized=true, vmin=-0.1, vmax=+0.1, cmap="seismic")
					clb = plt.colorbar(ima, ax=axs[2,:])
					clb.ax.set_ylabel("\$u^a-u^b\$")
					ima = axs[3,o].pcolormesh(transpose(au[z,:,:]), snap=true, shading="auto", rasterized=true, vmin=0.0, vmax=1.0, cmap="Oranges")
					clb = plt.colorbar(ima, ax=axs[3,:])
					clb.ax.set_ylabel("\$u^a\$")
				else
					axs[1,o].pcolormesh(transpose(gu[z,:,:]), snap=true, shading="auto", rasterized=true, vmin=0.0, vmax=1.0, cmap="Oranges")
					axs[2,o].pcolormesh(transpose(au[z,:,:]-gu[z,:,:]), snap=true, shading="auto", rasterized=true, vmin=-0.1, vmax=+0.1, cmap="seismic")
					axs[3,o].pcolormesh(transpose(au[z,:,:]), snap=true, shading="auto", rasterized=true, vmin=0.0, vmax=1.0, cmap="Oranges")
				end
			end
			plt.suptitle("\$t=$(t)\$ [ms]")
			plt.savefig(@sprintf("%s/%04d.png", frameDir, n), dpi=300, bbox_inches="tight")
		catch
			print("\n\t Stopping frame generation at iteration $(n).\n")
			break
		end
	end
	try
		rm("$(workDir)_innovation.mp4");
	catch
	end
	run(`ffmpeg -r 10 -f image2 -start_number 1 -i $(frameDir)/%04d.png -c:v h264_nvenc -pix_fmt yuv420p -loglevel quiet -stats $(workDir)_innovation.mp4`);
	rm(frameDir, recursive=true);
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
	#print("This is >7.23-1.25× faster than the Python implementation!\n")	
	innovationMovie(workDir)    
	print("\n")
end

main()
