ENV["TMPDIR"] = "/tmp"

using TorelliTrees
using Oscar

function same_poly_by_names(source_poly, target_poly)
    source_names = string.(gens(parent(source_poly)))
    target_names = string.(gens(parent(target_poly)))
    source_names == target_names || return false
    return evaluate(source_poly, gens(parent(target_poly))) == target_poly
end

function first_diff_line(a::String, b::String)
    al = split(a, "\n")
    bl = split(b, "\n")
    n = min(length(al), length(bl))
    for i in 1:n
        al[i] == bl[i] || return (i, al[i], bl[i])
    end
    if length(al) != length(bl)
        extra_a = length(al) >= n + 1 ? al[n + 1] : "<no line>"
        extra_b = length(bl) >= n + 1 ? bl[n + 1] : "<no line>"
        return (n + 1, extra_a, extra_b)
    end
    return nothing
end

function main()
    families = [[2, 1], [2, 2], [2, 3], [2, 4], [2, 5]]
    base = mktempdir()
    println("Working directory for comparison artifacts: ", base)
    println()

    all_ok = true
    for fam in families
        label = join(string.(fam), "_")
        prefix = "A" * join(string.(fam), "A")
        println("=== Family ", fam, " ===")

        opt_dir = joinpath(base, "opt_" * label)
        gen_dir = joinpath(base, "gen_" * label)
        mkpath(opt_dir)
        mkpath(gen_dir)

        println("Building trees for optimized run...")
        STs_opt = TorelliTrees.stratum_trees(fam)
        println("Tree count: ", length(STs_opt))
        println("Running compute_contributions_backup (optimized/default path)...")
        TorelliTrees.compute_contributions_backup(STs_opt, joinpath(opt_dir, "backup"))
        println("Generating optimized admcycles hodge code...")
        TorelliTrees.adm_code_hodge(STs_opt, prefix, joinpath(opt_dir, prefix))
        opt_code = read(joinpath(opt_dir, prefix * ".sage"), String)

        println("Building trees for forced-general run...")
        STs_gen = TorelliTrees.stratum_trees(fam)
        stage_ok = true
        mismatch_stage = nothing
        for (j, st) in enumerate(STs_gen)
            cont = TorelliTrees._compute_contribution(st; low_degree_optimization=false)
            opt_cont = STs_opt[j].cont_t[1].cont_t
            if !same_poly_by_names(cont, opt_cont)
                stage_ok = false
                mismatch_stage = (j, string(st), string(cont), string(opt_cont))
                break
            end
            if j % 10 == 0 || j == length(STs_gen)
                println("  Compared cont_t stage ", j, " / ", length(STs_gen))
            end
        end

        if !stage_ok
            all_ok = false
            println("STAGE COMPARISON FAILED for family ", fam)
            println("First mismatch at tree index ", mismatch_stage[1])
            println("Tree: ", mismatch_stage[2])
            println("Forced-general cont_t: ", mismatch_stage[3])
            println("Optimized/default cont_t: ", mismatch_stage[4])
            println()
            continue
        end

        println("Generating forced-general admcycles hodge code...")
        TorelliTrees.adm_code_hodge(STs_gen, prefix, joinpath(gen_dir, prefix))
        gen_code = read(joinpath(gen_dir, prefix * ".sage"), String)

        code_ok = opt_code == gen_code
        if !code_ok
            all_ok = false
            diff = first_diff_line(opt_code, gen_code)
            println("ADM CODE COMPARISON FAILED for family ", fam)
            if diff !== nothing
                println("First differing line: ", diff[1])
                println("Optimized/default: ", diff[2])
                println("Forced-general:   ", diff[3])
            end
        else
            println("Family ", fam, ": cont_t agrees at every stage and adm_code_hodge output matches exactly.")
        end
        println()
    end

    println("=== Final summary ===")
    println(all_ok ? "All requested comparisons passed." : "At least one comparison failed.")
    println("Artifacts kept in: ", base)
end

main()
