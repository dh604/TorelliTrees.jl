# Return (a,b) where a is the number of all g extremal trees and b is the number of
# irreducible ones.
function count_stratum_trees(g::Vector{Int64}; printStuff::Bool=false)::Tuple{Int64, Int64}
  ST = stratum_trees(g; printStuff=printStuff, compute_smoothings=false)
  n_st = length(ST)
  n_irred = 0
  for st in ST
    any(i -> i == 0, st.gen) && continue
    n_irred += 1
  end
  return (n_irred, n_st)
end

# Save the number of irreducible and any g extremal trees to the given file, for each of the entries g of gs.
function save_stratum_trees(gs::Vector{Vector{Int64}}, filename::String; printStuff::Bool=false)
  for g in gs
    a, b = count_stratum_trees(g; printStuff=printStuff)
    f = open(filename, "a")
    write(f, "Number of (irred/any) g-extremal trees for g=$g: $((a, b))\n")
    close(f)
    GC.gc()
  end
end