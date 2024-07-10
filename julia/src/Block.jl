# This file contains code that helps to query from the registry to determine the format

## Code dealing with Blocks 
@with_kw struct Block
    IMAX::Int64
    JMAX::Int64
    KMAX::Int64
    X::Array{Float64,3}
    Y::Array{Float64,3}
    Z::Array{Float64,3}
end

""" Reads a chunk X, Y, or Z variable from a plot3d file
    Args
        f: IO
        IMAX: 
        JMAX:
        KMAX
    Returns
        Array 3 dimensions for either X,Y, or Z
    """
function read_plot3D_chunk(f::IO,IMAX::Int64,JMAX::Int64,KMAX::Int64)::Array{Float64,3}
    A = Array{Float64,3}(undef,IMAX,JMAX,KMAX)
    for k in 1:KMAX
        for j in 1:JMAX
            for i in 1:IMAX
                A[i,j,k] = Float64(read(f,Float32))
            end
        end
    end
    return A
end


""" Reads a multi-block Plot3D File
    Args 
        filename (string)
"""
function read_multi_block(filename::String)
    
    blocks = Block[]
    if isfile(filename)
        io = open(filename,"r")
        IMAX = Int64[]; JMAX = Int64[]; KMAX = Int64[]
        NBLOCKS = Int64(read(io,Int32))
        for i in 1:NBLOCKS
            append!(IMAX,Int64(read(io,Int32)))
            append!(JMAX,Int64(read(io,Int32)))
            append!(KMAX,Int64(read(io,Int32)))
        end
        for i in 1:NBLOCKS
            X = read_plot3D_chunk(io,IMAX[i],JMAX[i],KMAX[i])
            Y = read_plot3D_chunk(io,IMAX[i],JMAX[i],KMAX[i])
            Z = read_plot3D_chunk(io,IMAX[i],JMAX[i],KMAX[i])
            b = Block(IMAX=IMAX[i],JMAX=JMAX[i],KMAX=KMAX[i],X=X,Y=Y,Z=Z)
            push!(blocks,b)
        end
        close(io)
    end
    return blocks
end
## End Blocks 