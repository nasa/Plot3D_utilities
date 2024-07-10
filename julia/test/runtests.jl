using Plot3D
using Test

@testset "test_read.jl" begin
    include("test_read.jl")
    @test read_VSPT() == "completed"
end
