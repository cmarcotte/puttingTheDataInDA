module readData

export readRestarts, readRestarts!, nx, ny, nz

const nx,ny,nz,nnz = 50,200,200,50
const u_tmp = zeros(Float64, (nnz+2), (ny+2), (nx+2));
const v_tmp = zeros(Float64, (nnz+1), (ny+1), (nx+1));
const w_tmp = zeros(Float64, (nnz+1), (ny+1), (nx+1));

function readRestarts(fnames; nproc=4)
	# this one just allocates, calls the in-place, and returns
	u = zeros(Float64,nx,ny,nz)
	v = zeros(Float64,nx,ny,nz)
	w = zeros(Float64,nx,ny,nz)
	
	readRestarts!(u, v, w, fnames; nproc=nproc)
	
	return (u,v,w)

end
function readRestarts!(u, v, w, fnames; nproc=4)
	# Check if the files and MPI processes match up first
	@assert length(fnames) == nproc
	# Check if the arrays are the correct size
	@assert size(u) == (nx,ny,nz)
	@assert size(v) == (nx,ny,nz)
	@assert size(w) == (nx,ny,nz)
	
	for (p,fname) in enumerate(fnames)
		open(fname,"r") do io
			gb = read(io, Float32)
			read!(io, u_tmp)
			read!(io, v_tmp)
			read!(io, w_tmp)
			u[:,:,1+(p-1)*nnz:(p+0)*nnz] .= u_tmp[begin+1:end-1,begin+1:end-1,begin+1:end-1]
			v[:,:,1+(p-1)*nnz:(p+0)*nnz] .= v_tmp[begin+1:end-0,begin+1:end-0,begin+1:end-0]
			w[:,:,1+(p-1)*nnz:(p+0)*nnz] .= w_tmp[begin+1:end-0,begin+1:end-0,begin+1:end-0]
		end
	end
	return nothing

end

function readRestarts(fnames, index::Array{T}=[1]; nproc=4) where T <: Int 
	# honestly this is probably not useful. Easier to just discard unneeded arrays that are read in anyway.
	# this one needs to be allocating because the number of fields is known at runtime
	# this one only loads a subset of the fields
	@assert index == unique(index)
	arrs = [zeros(Float64,nx,ny,nz) for i in index]
	
	for (p,fname) in enumerate(fnames)
		open(fname,"r") do io
			gb = read(io, Float32)
			read!(io, u_tmp)
			read!(io, v_tmp)
			read!(io, w_tmp)
			if 1 in index
				arrind = findfirst(index.==1)
				arrs[arrind][:,:,1+(p-1)*nnz:(p+0)*nnz] .= u_tmp[begin+1:end-1,begin+1:end-1,begin+1:end-1]
			elseif 2 in index
				arrind = findfirst(index.==2)
				arrs[arrind][:,:,1+(p-1)*nnz:(p+0)*nnz] .= v_tmp[begin+1:end-0,begin+1:end-0,begin+1:end-0]
			elseif 3 in index
				arrind = findfirst(index.==3)
				arrs[arrind][:,:,1+(p-1)*nnz:(p+0)*nnz] .= w_tmp[begin+1:end-0,begin+1:end-0,begin+1:end-0]
			end
		end
	end
	return arrs
end

end

using Main.readData, WGLMakie, Printf
# Original testing of Makie:

u = rand(Float64,nx,ny,nz)
v = rand(Float64,nx,ny,nz)
w = rand(Float64,nx,ny,nz)
U = Observable(u)
V = Observable(v)
W = Observable(w)

fig = Figure(resolution = (1200, 800))
axs = [Axis3(fig[1,i], aspect=(0.25,1,1)) for i in 1:3]
volume!(axs[1],U,alpha=0.02,colormap="Oranges")
volume!(axs[2],V,alpha=0.02,colormap="Purples")
volume!(axs[3],W,alpha=0.02,colormap="Greens")
for ax in axs
	hidedecorations!(ax)
end
Colorbar(fig[2,1]; vertical = false, flipaxis = false, colorrange=(0,1), colormap="Oranges", label=L"$u$")
Colorbar(fig[2,2]; vertical = false, flipaxis = false, colorrange=(0,1), colormap="Purples", label=L"$v$")
Colorbar(fig[2,3]; vertical = false, flipaxis = false, colorrange=(0,1), colormap="Greens",  label=L"$w$")
tit = Label(fig[0, :], text = L"$t = 0$ [ms]", textsize = 32)
resize_to_layout!(fig)

#contour(fig[1,1],U,alpha=0.05,colormap="Oranges")
#contour(fig[1,2],V,alpha=0.05,colormap="Purples")
#contour(fig[1,3],W,alpha=0.05,colormap="Greens")

#U = Observable(u[1,:,:])
#V = Observable(v[1,:,:])
#W = Observable(w[1,:,:])
#heatmap(fig[1,1],U,colormap="Oranges")
#heatmap(fig[1,2],V,colormap="Purples")
#heatmap(fig[1,3],W,colormap="Greens")

#volumeslices(fig[1,1],1:50,1:200,1:200,U,colormap="Oranges",colorrange=(0,1),transparency=true)
#volumeslices(fig[1,2],1:50,1:200,1:200,V,colormap="Purples",colorrange=(0,1),transparency=true)
#volumeslices(fig[1,3],1:50,1:200,1:200,W,colormap="Greens",colorrange=(0,1),transparency=true)

record(fig, "volume.mp4", 0:124; framerate=20) do t
	fnames = [@sprintf("./%04d/restart3d.%03d",t,p) for p in 0:3]
	readRestarts!(u, v, w, fnames)
	tit.text[] = L"$t = %$(5.0*t)$ [ms]"
	U[] = u#[1,:,:]
	V[] = v#[1,:,:]
	W[] = w#[1,:,:]
end


# Specifically heatmap testing of CairoMakie to get some 3D info while headless
using CairoMakie, Printf

u = zeros(Float64,nx,ny,nz)
v = zeros(Float64,nx,ny,nz)
w = zeros(Float64,nx,ny,nz)

Zs = [1,5,10,15,20,25,30,35,40,35,50]

Us = [Observable(u[z,:,:]) for z in Zs]
Vs = [Observable(v[z,:,:]) for z in Zs]
Ws = [Observable(w[z,:,:]) for z in Zs]

fig = Figure(resolution = (200*length(Zs), 200*3))
axs = [ Axis(fig[i, j], width = ny, height = nz) for i in 1:3, j in 1:length(Zs)]

for ax in axs
	hidedecorations!(ax)
end

for j in 1:length(Zs)
	heatmap!(axs[1,j],Us[j],colormap="Oranges",colorrange=(0,1))
	heatmap!(axs[2,j],Vs[j],colormap="Purples",colorrange=(0,1))
	heatmap!(axs[3,j],Ws[j],colormap="Greens", colorrange=(0,1))
	axs[1,j].title = "z=$(Zs[j])"
end
Colorbar(fig[1,length(Zs)+1]; colorrange=(0,1), colormap="Oranges", label=L"u")
Colorbar(fig[2,length(Zs)+1]; colorrange=(0,1), colormap="Purples", label=L"v")
Colorbar(fig[3,length(Zs)+1]; colorrange=(0,1), colormap="Greens",  label=L"w")

resize_to_layout!(fig)
save("heatmap_array.png", fig)

record(fig, "heatmaparray.mp4", 0:124; framerate=20) do t
	fnames = [@sprintf("./%04d/restart3d.%03d",t,p) for p in 0:3]
	readRestarts!(u, v, w, fnames)
	for (j,z) in enumerate(Zs)
		Us[j][] = u[z,:,:]
		Vs[j][] = v[z,:,:]
		Ws[j][] = w[z,:,:]
	end
end

end
