module Plot3D    
    using Parameters, IterTools, Printf
    using ProgressMeter
    using DataFrames
    import Base

    include("Block.jl")
    include("Face.jl")
    export Block, read_multi_block              # From Block
    export Face, match_indices, add_vertex      # From Face
    export find_block_connectivity


    ## Internal Functions related to find_block_connectivity
    function get_block_outer_faces(block::Block)::Array{Face,1}
        """ Get the outer faces of a block

                Args
                    block1: Block information
                Returns
                    List of outer faces (Face)
            """
        I = [1,block.IMAX]
        J = [1,block.JMAX]
        K = [1,block.KMAX]
            # Create outer faces
        faces = Face[]

        face = Face()
        i = I[1]
        for j in J
            for k in K
                add_face_vertex(face, block.X[i,j,k], block.Y[i,j,k], block.Z[i,j,k], i, j, k)
            end
        end
        push!(faces, face)

        face = Face()
        i = I[2]
        for j in J
            for k in K
                add_face_vertex(face, block.X[i,j,k], block.Y[i,j,k], block.Z[i,j,k], i, j, k)
            end
        end
        push!(faces, face)

        face = Face()
        j = J[1]
        for i in I
            for k in K
                add_face_vertex(face, block.X[i,j,k], block.Y[i,j,k], block.Z[i,j,k], i, j, k)
            end
        end
        push!(faces, face)

        face = Face()
        j = J[2]
        for i in I
            for k in K
                add_face_vertex(face, block.X[i,j,k], block.Y[i,j,k], block.Z[i,j,k], i, j, k)
            end
        end
        push!(faces, face)

        face = Face()
        k = K[1]
        for i in I
            for j in J
                add_face_vertex(face, block.X[i,j,k], block.Y[i,j,k], block.Z[i,j,k], i, j, k)
            end
        end
        push!(faces, face)

        face = Face()
        k = K[2]
        for i in I
            for j in J
                add_face_vertex(face, block.X[i,j,k], block.Y[i,j,k], block.Z[i,j,k], i, j, k)
            end
        end
        push!(faces, face)

            # Check if faces match each other
        non_matching = Face[]
        for i in 1:length(faces)
            matching = true
            for j in 1:length(faces)
                if (i != j) && (faces[i] != faces[j])
                    matching = false
                end
            end
            if !matching
                push!(non_matching, faces[i]) # these are guaranteed to be exterior
            end
        end
        return non_matching # these should be the outer faces
    end

    function point_match(x::Float64, y::Float64, z::Float64, x2::Array{Float64,3}, y2::Array{Float64,3}, z2::Array{Float64,3})
        """ Checks to see if x,y,z is present in a block
            Args
                x coordinates to check if inside block
                y
                z
                x2 all the x coordinates of block 2
                y2 all the y coordinates of block 2
                z2 all the z coordinates of block 2
            Return
                True/False: checks if there's a match
                i,j,k: index in block where match is found
            """
        dx = x .- x2
        dy = y .- y2
        dz = z .- z2
        d = sqrt.(dx.* dx .+ dy.* dy .+ dz.* dz)
        val, location = findmin(d)

        if val < 1E-6
            return true, location[1], location[2], location[3]
        end
        return false, 0, 0, 0
    end

    function get_face_intersection(face1::Face, face2::Face, block1::Block, block2::Block)
        """Get the index of the intersection between two faces located on two different blocks

                Args
                    face1: A face, exterior face hopefully
                    face2: A face, exterior face hopefully
                    block1: block for face1
                    block2: block for face2

                Returns
                    match_locations
                        A list of all indices in block1 that match block2
                            [ [(block1_i,block1_j,block1_k), (block2_i, block2_j, block2_k)],
                                [(block1_i,block1_j,block1_k), (block2_i, block2_j, block2_k)] ]

                    match_location_corners
                        A dictionary with the intersections
                        Block1: {I: [IMIN,IMAX], J: [JMIN,JMAX], K: [KMIN,KMAX]}
                        Block2: {I: [IMIN,IMAX], J: [JMIN,JMAX], K: [KMIN,KMAX]}
            """
        I1 = [minimum(face1.I),maximum(face1.I)]
        J1 = [minimum(face1.J),maximum(face1.J)]
        K1 = [minimum(face1.K),maximum(face1.K)]
        # @printf("Block1 Range = %3d:%3d %3d:%3d %3d:%3d\n",I1[1],I1[2],J1[1],J1[2],K1[1],K1[2])

        I2 = [minimum(face2.I),maximum(face2.I)]
        J2 = [minimum(face2.J),maximum(face2.J)]
        K2 = [minimum(face2.K),maximum(face2.K)]
        # @printf("Block2 Range = %3d:%3d %3d:%3d %3d:%3d\n",I2[1],I2[2],J2[1],J2[2],K2[1],K2[2])

        x2 = block2.X[I2[1]:I2[2],J2[1]:J2[2],K2[1]:K2[2]]
        y2 = block2.Y[I2[1]:I2[2],J2[1]:J2[2],K2[1]:K2[2]]
        z2 = block2.Z[I2[1]:I2[2],J2[1]:J2[2],K2[1]:K2[2]]
        """
            General Search
        """
        match_locations = Any[]        
        if I1[1] == I1[2]
            i = I1[1]
            combo = IterTools.product(J1[1]:J1[2], K1[1]:K1[2])
            for p in combo
                j, k = p
                x = block1.X[i,j,k]
                y = block1.Y[i,j,k]
                z = block1.Z[i,j,k]
                bIntersect, i2, j2, k2 = point_match(x, y, z, x2, y2, z2)
                if bIntersect
                    push!(match_locations, [(i, j, k),(i2, j2, k2)]) 
                end
            end
        end

        if J1[1] == J1[2]
            j = J1[1]
            combo = IterTools.product(I1[1]:I1[2],K1[1]:K1[2])
            for p in combo
                i, k = p
                x = block1.X[i,j,k]
                y = block1.Y[i,j,k]
                z = block1.Z[i,j,k]
                bIntersect, i2, j2, k2 = point_match(x, y, z, x2, y2, z2)
                if bIntersect
                    push!(match_locations, [(i, j, k),(i2, j2, k2)]) 
                end
            end
        end

        if K1[1] == K1[2]
            k = K1[1]
            combo = IterTools.product(I1[1]:I1[2], J1[1]:J1[2])
            for p in combo
                i, j = p
                x = block1.X[i,j,k]
                y = block1.Y[i,j,k]
                z = block1.Z[i,j,k]
                bIntersect, i2, j2, k2 = point_match(x, y, z, x2, y2, z2)
                if bIntersect
                    push!(match_locations, [(i, j, k),(i2, j2, k2)]) 
                end
            end
        end

        match_location_corners = Any[]
        # Find matched the corners
        df_block1 = DataFrame(i = match_locations[:][1][1], j =  match_locations[:][1][2], k =  match_locations[:][1][3]);
        df_block2 = DataFrame(i = match_locations[:][2][1], j =  match_locations[:][2][2], k =  match_locations[:][2][3]);

        if I1[1] == I1[2]
            sort(df_block1, [:j :k])
        elseif J1[1] == J1[2]
            sort(df_block1, [:i :k])
        elseif K1[1] == K1[2]               # Bottom Corner IMIN,JMIN Top Corner IMAX,JMAX
            sort(df_block1, [:i :j])
        end
        
        if I2[1] == I2[2]
            sort(df_block2, [:j :k])
        elseif J1[1] == J1[2]
            sort(df_block2, [:i :k])
        elseif K1[1] == K1[2]
            sort(df_block2, [:i :j])
        end
        
        b1_corners = Dict("block1" => Dict("IMIN" => I1[1], "JMIN" => J1[1], "KMIN" => K1[1],
                                        "IMAX" => I1[2], "JMAX" => J1[2], "KMAX" => K1[2],
                                        "corner1" => vec(df_block1[1,:]), "corner2" => vec(df_block1[end,:])))

        b2_corners = Dict("block2" => Dict("IMIN" => I2[1], "JMIN" => J2[1], "KMIN" => K2[1],
                                        "IMAX" => I2[2], "JMAX" => J2[2], "KMAX" => K2[2],
                                        "corner1" => vec(df_block2[end,:]), "corner2" => vec(df_block2[end,:])))
        
                                                
        

        return match_locations, (b1_corners, b2_corners)
    end

    function find_block_connectivity(block1::Block, block2::Block)
        """Takes two blocks and finds what indices match 
                Looks at exterior faces of block1 and block2.
                Checks to see if any of the exterior faces match each other 

                Args:
                    block1: a plot3d block
                    block2: a plot3d block

                Returns:
                    corners of matching pair as block1_corners,block2_corners
                    ([imin,jmin,kmin],[imax,jmax,kmax]), ([imin,jmin,kmin],[imax,jmax,kmax])
        """
        # Check to see if outer face of block 1 matches any of the outer faces of block 2
        block_match_indices = Any[]
        block_match_corners = Any[]

        block1_outer = get_block_outer_faces(block1)
        block2_outer = get_block_outer_faces(block2)
            # Checks the nodes of the outer faces to see if any of them match

        print("\nLooping through all exterior faces to find a match\n")
        @printf("Block 1 Face | Block 2 Face\n")
        match_connections, matched_corners = get_face_intersection(block1_outer[4], block2_outer[5], block1, block2)
        push!(block_match_indices, match_connections)
        push!(block_match_corners, matched_corners)

        # for p in 1:length(block1_outer)
        #     block1_face = block1_outer[p]            
        #     for q in 1:length(block2_outer)
        #         @printf("     %2d     |      %2d    \n",p,q)
        #         block2_face = block2_outer[q]
        #         match_connections, matched_corners = get_face_intersection(block1_face,block2_face,block1,block2)
        #         if length(match_connections)>0
        #             push!(block_match_indices, match_connections)
        #             push!(block_match_corners, matched_corners)
        #         end
        #     end
        # end
        return block_match_indices, block_match_corners
    end


end