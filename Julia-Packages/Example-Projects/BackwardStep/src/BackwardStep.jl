using Plot3D: read_multi_block, Block, find_block_connectivity
using JLD


## Cell 1 Reads the block data 
if !isfile("blocks.jld")
    blocks = read_multi_block("../../../../testfiles/8block.p3d")
    jldopen("blocks.jld", "w") do file
        write(file, "blocks", blocks)  # alternatively, say "@write file A"
    end
end
## 
blocks = jldopen("blocks.jld", "r") do file
    read(file, "blocks")
end

## Cell 2 Read the Block data 
if !isfile("block_matches.jld")
    print("Block read, Finding matching faces\n")
    block_match_indices, block_match_corners = find_block_connectivity(blocks[1],blocks[2])
    jldopen("block_matches.jld", "w") do file        
        write(file, "block_match_indices", block_match_indices)  # alternatively, say "@write file A"
        write(file, "block_match_corners", block_match_corners)  # alternatively, say "@write file A"
    end
end

block_matches = jldopen("block_matches.jld", "r") do file
    read(file, "block_match_indices")
end
print("check")
##
