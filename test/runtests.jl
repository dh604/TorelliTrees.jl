using TorelliTrees
using Test

@testset "TorelliTrees.jl" begin
    # Test the number of graphs that we get
    T13 = stratum_trees([1, 3])
    T14 = stratum_trees([1, 4])
    T15 = stratum_trees([1, 5])
    T24 = stratum_trees([2, 4])
    T33 = stratum_trees([3, 3])
    @test length(T13) == 4
    @test length(T14) == 10
    @test length(T15) == 24
    @test length(T24) == 153
    @test length(T33) == 210
    @test count(st -> length(st.min_smoothings) == 0, T13) == 3
    @test count(st -> length(st.min_smoothings) == 0, T14) == 5
    @test count(st -> length(st.min_smoothings) == 0, T15) == 7
    @test count(st -> length(st.min_smoothings) == 0, T24) == 15
    @test count(st -> length(st.min_smoothings) == 0, T33) == 19
    
    # Test if calculating contributions works
    compute_contributions(T13)
    @test true
    compute_contributions(T14)
    @test true
    compute_contributions(T15)
    @test true
    compute_contributions(T24)
    @test true
    compute_contributions(T33)
    @test true

    # Test if writing the admcycles code works
    adm_code(T13, "A1A3", "testoutput/A1A3")
    @test true
    adm_code(T14, "A1A4", "testoutput/A1A4")
    @test true
    adm_code(T15, "A1A5", "testoutput/A1A5")
    @test true
    adm_code(T24, "A2A4", "testoutput/A2A4")
    @test true
    adm_code(T33, "A3A3", "testoutput/A3A3")
    @test true
end
