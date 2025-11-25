function compute_contributions(STs::Vector{StratumTree}; print_contributions::Bool=false)
  for st in STs
    cont = _compute_contribution(st)
    !print_contributions && continue
    println("Contribution of $st:")
    println(cont)
    println()
  end
end

# Compute the contribution polynomial of this stratum tree, assuming the contributions of all its smoothings have been computed.
function _compute_contribution(st::StratumTree)
  ne = Graphs.ne(st.GG)
  cps = _critical_paths(st)
  N_deg = _tor_pullback_deg(st) # >= _cont_deg(st), so there are always enough chern classes to substitute for the l's.
  nl = N_deg - length(cps)
  cont_t_deg = _cont_deg(st)
  
  # Create ring in which Cont_T lives
  edge_vect = collect(Graphs.edges(st.GG))
  edge_to_gen = Dict{GraphsEdgeType, PolyType}()
  var_names = vcat(["x$(Graphs.src(e))$(Graphs.dst(e))" for e in edge_vect], ["ls$i" for i in 1:nl], ["c$i" for i in 1:N_deg])
  var_degs = vcat([1 for _ in edge_vect], [i for i in 1:nl], [i for i in 1:N_deg])
  R, g = graded_polynomial_ring(QQ, var_names, var_degs)
  ls = g[ne+1:ne+nl]
  ci = g[ne+nl+1:ne+nl+N_deg]
  for j in 1:ne
    e = edge_vect[j]
    edge_to_gen[e] = g[j]
    edge_to_gen[Graphs.reverse(e)] = g[j]
  end
  c_tot_z_l = _tot_chern_class(edge_to_gen, ls, cps)
  c_tot_c = sum(vcat([one(R)], ci))
  # println(factor(c_tot_z_l))
  # println(c_tot_c)

  if count(c -> c == 1, st.col) == 0
    # Irreducible component, apply excess intersection formula
    denom_factors = [1 + edge_to_gen[e] for e in edge_vect]
    # result in terms of z_e's and ls_i's.
    cont_t = homogeneous_component(c_tot_c * _denominator_to_power_series(denom_factors, cont_t_deg), cont_t_deg)
    st.cont_t[1] = Cont_t_container(cont_t, ne, nl, N_deg, edge_vect, edge_to_gen)
    return cont_t
  else
    # Not an irreducible component, apply the recursion.
    cont_t_times_z = ci[N_deg]
    for sm in st.all_smoothings
      container = sm.ST1.cont_t[1]
      cont_of_sm = container.cont_t
      eval_vector = vcat([edge_to_gen[sm.edge_map[e]] for e in container.edge_vect], [zero(R) for _ in 1:container.nl], ci[1:container.N_deg])
      # println(eval_vector)
      cont_of_sm_evaluated = evaluate(cont_of_sm, eval_vector)
      cont_t_times_z -= cont_of_sm_evaluated * prod(e -> edge_to_gen[sm.edge_map[e]], container.edge_vect)
    end
    # println("Cont_t_times z: $cont_t_times_z")

    # First, we need to express everything in terms of ls.
    eval_vec_ci_to_ls = copy(gens(R))
    foreach(i -> eval_vec_ci_to_ls[ne+nl+i] = homogeneous_component(c_tot_z_l, i), 1:N_deg)
    cont_t_times_z = evaluate(cont_t_times_z, eval_vec_ci_to_ls)

    # Second, we need to substitute the (fewer) ls in terms of chern classes again.
    l_tot_z_c = c_tot_c * _denominator_to_power_series(_crit_path_factors(edge_to_gen, cps), nl) # was N_deg
    eval_vec_ls_to_ci = copy(gens(R))
    foreach(i -> eval_vec_ls_to_ci[ne+i] = homogeneous_component(l_tot_z_c, i), 1:nl)
    cont_t_times_z = evaluate(cont_t_times_z, eval_vec_ls_to_ci)

    # Finally, we need to divide the result by the edge product and are done.
    prod_e = prod(e -> edge_to_gen[e], edge_vect)
    @req is_divisible_by(cont_t_times_z, prod_e) "cont_t_times_z is not divisible by prod_e!"
    cont_t = div(cont_t_times_z, prod_e)
    st.cont_t[1] = Cont_t_container(cont_t, ne, nl, N_deg, edge_vect, edge_to_gen)
    return cont_t
  end
end

function _crit_path_prod(edge_to_gen::Dict{GraphsEdgeType, PolyType}, cps::Vector{Vector{GraphsEdgeType}})::PolyType
  # res = one(edge_to_gen[cps[1][1]])
  # # factor for each path
  # for p in cps
  #   tmp = one(res)
  #   for e in p
  #     tmp += edge_to_gen[e]
  #   end
  #   res *= tmp
  # end
  # return res
  return prod(_crit_path_factors(edge_to_gen, cps))
end

function _crit_path_factors(edge_to_gen::Dict{GraphsEdgeType, PolyType}, cps::Vector{Vector{GraphsEdgeType}})::Vector{PolyType}
  R = parent(edge_to_gen[cps[1][1]])
  return [sum(vcat([one(R)], [edge_to_gen[e] for e in p])) for p in cps]
end

function _tot_chern_class(edge_to_gen::Dict{GraphsEdgeType, PolyType}, ls::Vector{PolyType}, cps::Vector{Vector{GraphsEdgeType}})::PolyType
  res = _crit_path_prod(edge_to_gen, cps)
  # l factors
  tmp = one(res)
  for lsi in ls
    tmp += lsi
  end
  res *= tmp
  return res
end

# This assumes that f-1 is of homogeneous degree 1.
# When this assumption is not met, more term than required are returned (which is fine for our calculations).
function __denominator_to_power_series(f::PolyType, maxDeg::Int64)::PolyType
  x = 1-f
  res = one(f)
  fLast = one(f)
  for _ in 1:maxDeg
    fLast *= x
    res += fLast
  end
  return res
end

# Here each entry of f should be 1 + homogeneous part of degree one.
function _denominator_to_power_series(f::Vector{PolyType}, maxDeg::Int64)::PolyType
  res = __denominator_to_power_series(f[1], maxDeg)
  for i in 2:length(f)
    res *= __denominator_to_power_series(f[i], maxDeg)
    res = sum(j -> homogeneous_component(res, j), 0:maxDeg)
  end
  return res
end

# Use Oscar's homogeneous_component instead.
# function _homogeneous_part(f::PolyType, d::Int64)::QQMPolyRingElem
#   R = parent(f)
#   n = number_of_variables(R)
#   g = gens(R)
#   res = zero(R)
#   for e in Oscar.weak_compositions(d, n)
#     e = [i for i in e]
#     res += coeff(f, e) * prod(g .^ e)
#   end
#   return res
# end

# The degree of Cont_T, which is the degree of Tor^*(...) minus the number of edges.
function _cont_deg(st::StratumTree)::Int64
  return _tor_pullback_deg(st) - Graphs.ne(st.GG)
end

# The degree of the final pullback class Tor^*(...).
function _tor_pullback_deg(st::StratumTree)::Int64
  nv = Graphs.nv(st.GG)
  g_tot = sum(st.gen)
  n_cols = maximum(st.col)
  g = [sum([st.gen[v] for v in 1:nv if st.col[v] == c]) for c in 1:n_cols]
  return div(g_tot*(g_tot+1) - sum(gk -> gk*(gk+1), g), 2)
end

# Return the critical paths of the stratum tree (only one choice of direction per path.)
function _critical_paths(st::StratumTree)::Vector{Vector{GraphsEdgeType}}
  res = Vector{Vector{GraphsEdgeType}}()
  nv = Graphs.nv(st.GG)
  for v in 1:nv
    for w in v+1:nv
      (st.col[v] == st.col[w] || st.col[v] == 1 || st.col[w] == 1) && continue
      # this picks the shortest path, which is automatically the unique path.
      path = Graphs.a_star(st.GG, v, w)
      lp = length(path)
      # Require that the path has only genus zero interior vertices.
      lp > 1 && any(j -> st.col[Graphs.src(path[j])] != 1, 2:lp) && continue
      # The path is critical as it has only genus zero interior vertices and connects two positive genus vertices of different colors.
      push!(res, path)
    end
  end
  return res
end