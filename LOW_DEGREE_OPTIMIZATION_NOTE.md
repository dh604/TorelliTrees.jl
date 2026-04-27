# Low-degree contribution optimization

This note summarizes the optimization added for contribution computations in the family `[2, g-2]`, together with the current reproducibility script.

## Idea

For a stratum tree `T`, the contribution polynomial `Cont_T` has degree

`deg(Cont_T) = N - e`,

where `N` is the Torelli pullback degree and `e` is the number of edges.

For the family `[2, g-2]`, this becomes

`N = 2g - 4`, so `deg(Cont_T) = 2g - 4 - e`.

Hence:

- `Cont_T` is linear exactly when `e = 2g - 5`.
- `Cont_T` is constant exactly when `e = 2g - 4`.

In those cases, the full recursive computation carries more Chern-class variables than are actually needed in the final answer. The optimization therefore uses a reduced graded polynomial ring containing:

- the edge variables,
- all `ls_i`,
- only `c1`.

The recursion is then computed in three steps:

1. Substitute `c_i -> homogeneous_component(c_tot_z_l, i)` exactly, as in the general code.
2. Perform the low-degree version of the second substitution using `1 + c1` instead of the full total Chern class `1 + c1 + c2 + ...`.
3. Divide by the edge product and lift the result back into the usual full ring.

The reason this works is that in degree `0` or `1`, terms involving `c_i` for `i >= 2` cannot survive in the final contribution.

## Current activation rule

The low-degree branch is implemented in [src/contributions.jl](/home/daniel/julia_workspace/TorelliTrees.jl/src/contributions.jl) inside `_compute_contribution`.

At the moment, the code uses the optimization whenever

`cont_t_deg <= 1`.

There used to be an additional guard `&& N_deg >= 8`, but it is currently commented out to allow broader testing.

## Reproducibility script

A standalone comparison script now lives at [test/compare_low_degree.jl](/home/daniel/julia_workspace/TorelliTrees.jl/test/compare_low_degree.jl).

It compares the default path against the forced general path `low_degree_optimization=false` on the families

- `[2,1]`
- `[2,2]`
- `[2,3]`
- `[2,4]`
- `[2,5]`

For each family it checks:

1. The stage-by-stage `cont_t` values agree, up to identification of generators by name rather than literal ring identity.
2. The final `adm_code_hodge` output agrees exactly.

You can rerun it from the project root with:

```bash
julia --startup-file=no --project=. test/compare_low_degree.jl
```

## Verification summary

The comparison script was run on 2026-04-27 and all requested comparisons passed.

Tree counts by family:

- `[2,1]`: 2 trees
- `[2,2]`: 9 trees
- `[2,3]`: 37 trees
- `[2,4]`: 153 trees
- `[2,5]`: 622 trees

Outcome:

- The computed `cont_t` agreed at every stage for all five families.
- The generated `adm_code_hodge` code agreed exactly for all five families.
