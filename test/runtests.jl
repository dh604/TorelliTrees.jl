using TorelliTrees
using Graphs
using Oscar
using Test

@testset "TorelliTrees.jl" begin
    testoutput_dir = joinpath(@__DIR__, "testoutput")
    mkpath(testoutput_dir)

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

    # Test the low-degree optimization on the [2, g-2] family for g up to 6.
    function _same_poly_in_parent(source_poly, target_poly)
        return evaluate(source_poly, gens(parent(target_poly))) == target_poly
    end

    for fam in ([2, 1], [2, 2], [2, 3], [2, 4])
        STs = fam == [2, 4] ? T24 : stratum_trees(fam)
        fam == [2, 4] || compute_contributions(STs)
        total_genus = sum(fam)
        @test all(st -> TorelliTrees._cont_deg(st) == 2 * total_genus - 4 - Graphs.ne(st.GG), STs)
        for st in STs
            TorelliTrees._cont_deg(st) <= 1 || continue
            isempty(st.all_smoothings) && continue
            optimized = st.cont_t[1].cont_t
            general = TorelliTrees._compute_contribution(st; low_degree_optimization=false)
            @test _same_poly_in_parent(general, optimized)
        end
    end

    # Test if writing the admcycles code works
    adm_code(T13, "A1A3", joinpath(testoutput_dir, "A1A3"))
    @test true
    adm_code(T14, "A1A4", joinpath(testoutput_dir, "A1A4"))
    @test true
    adm_code(T15, "A1A5", joinpath(testoutput_dir, "A1A5"))
    @test true
    adm_code(T24, "A2A4", joinpath(testoutput_dir, "A2A4"))
    @test true
    adm_code(T33, "A3A3", joinpath(testoutput_dir, "A3A3"))
    @test true
end
