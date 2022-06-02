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

