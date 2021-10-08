## Code dealing with Face 
mutable struct Face
    nvertex::Int64
    X::Array{Float64,1}
    Y::Array{Float64,1}
    Z::Array{Float64,1}
    I::Array{Int64,1}
    J::Array{Int64,1}
    K::Array{Int64,1}  
end

"""Default Constructor for Face 
"""
function Face()    
    return Face(0,Float64[], Float64[],Float64[],Int64[],Int64[],Int64[])
end

""" Add a vertex to a face 
    Args
        f: Face 
        x: x-vertex
        y: y-vertex
        z: z-vertex
        i: index of x
        j: index of y
        k: index of z
"""
function add_face_vertex(f::Face, x::Float64, y::Float64, z::Float64,i::Int64,j::Int64,k::Int64)
    push!(f.X,x)
    push!(f.Y,y)
    push!(f.Z,z)
    push!(f.I,i)
    push!(f.J,j)
    push!(f.K,k)
    f.nvertex+=1
end


""" Check to see if two faces are the same by matching vertices
    Args
        f: another face
    Returns
        match_indices e.g. [[1,2],[1,4]] Face 1 matches Face 2, Face 1 matches Face 4
"""
function match_face_indicies(f1::Face,f2::Face)
    
    tol = 1E-6
    matchedIndices = Any[]
    for i in 1:f1.nvertex
       for j in 1:f1.nvertex
           dx = abs(f1.X[i] - f2.X[j])
           dy = abs(f1.Y[i] - f2.Y[j])
           dz = abs(f1.Z[i] - f2.Z[j])
           if ((dx<tol) && (dy<tol) && (dz<tol))
               push!(matchedIndices,[i j])
           end
        end
    end
    return matchedIndices
end

Base.:(==)(lhs::Face, rhs::Face) = (length(match_face_indicies(lhs,rhs)) == rhs.nvertex)
Base.:(!=)(lhs::Face, rhs::Face) = (length(match_face_indicies(lhs,rhs)) != lhs.nvertex)
## End Face