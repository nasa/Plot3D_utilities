using Plot3D: read_multi_block, Block, find_block_connectivity
using JLD2
using JSON

## Cell 1 Reads the block data 
if !isfile("blocks.jld")
    blocks = read_multi_block("../../../testfiles/finalmesh.xyz")
    @save "blocks.jld" blocks
end
## 
@load "blocks.jld" blocks

## Cell 2 Read the Block data 
print("Block read, Finding matching faces\n")
block_match_indices, block_match_corners = find_block_connectivity(blocks[1],blocks[2])
json_string = JSON.json(block_match_indices)
open("block_match_indices.json","w") do f
    JSON.print(f, json_string, 4)
end

json_string = JSON.json(block_match_corners)
open("block_match_corners.json","w") do f
    JSON.print(f, json_string)
end


# @load "block_match_indices.jld" block_match_indices
# @load "block_match_corners.jld" block_match_corners


print("check")
##

