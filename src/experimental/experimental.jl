#######################################################################################
# Some experimental code to generalize and test Pixton's formula.                     #
# The code in this file is not used outside of the "experimental" folder.             #
#######################################################################################

function _set_chern_to_one(cont::Cont_t_container)
  R = parent(cont.cont_t)
  subst_vect = vcat(gens(R)[1:cont.ne+cont.nl], [one(R) for _ in 1:cont.N_deg])
  cont_subst = evaluate(cont.cont_t, subst_vect)
  return cont_subst
end

function _pixton_formula(st::StratumTree)::PolyType
  cont = st.cont_t[1]
  R = parent(cont.cont_t)
  res = one(R)
  nv = Graphs.nv(st.GG)
  # println("nv=$nv")
  root = 0
  for v in 1:nv
    if st.col[v] == 2
      root = v
      break
    end
  end
  # println("Root: $root")
  # @req root != 0 "No root was found"
  denoms = Vector{PolyType}()
  for v in 1:nv
    # println("vert $v has col $(st.col[v]) and gen $(st.gen[v])")
    if st.gen[v] == 0
      p = Graphs.a_star(st.GG, v, root)
      res *= (sum(vcat([one(R)], [cont.edge_to_gen[e] for e in p])))^(Graphs.degree(st.GG, v)-2)
    elseif st.col[v] == 3
      p = Graphs.a_star(st.GG, v, root)
      path_sum = sum(vcat([one(R)], [cont.edge_to_gen[e] for e in p]))
      push!(denoms, path_sum)
    end
  end
  # println("Denoms: $denoms")
  res *= _denominator_to_power_series(denoms, cont.N_deg)
  res = div(res, prod([cont.edge_to_gen[e] for e in cont.edge_vect]))
  res = sum([homogeneous_component(res, i) for i in 0:cont.N_deg-cont.ne])
  res = res * (-1)^count(v -> st.col[v] == 3, 1:nv)
  return res
end

function _pixton_generalization_1(st::StratumTree)::PolyType
  cont = st.cont_t[1]
  R = parent(cont.cont_t)
  res = one(R)
  nv = Graphs.nv(st.GG)
  denoms = Vector{PolyType}()
  for p in _critical_paths(st)
    path_sum = sum(vcat([one(R)], [cont.edge_to_gen[e] for e in p]))
    push!(denoms, path_sum)
  end
  for v in 1:nv
    # println("vert $v has col $(st.col[v]) and gen $(st.gen[v])")
    st.gen[v] != 0 && continue
    for w in 1:nv
      # st.col[w] != 2 && continue
      st.col[w] != 2 && continue
      p = Graphs.a_star(st.GG, v, w)
      any(e -> st.col[Graphs.src(e)] != 1, p) && continue
      res *= (sum(vcat([one(R)], [cont.edge_to_gen[e] for e in p])))^(Graphs.degree(st.GG, v)-2)
      # res = sum([homogeneous_component(res, i) for i in 0:cont.N_deg-cont.ne])
    end
  end
  # println("Denoms: $denoms")
  res *= _denominator_to_power_series(denoms, cont.N_deg)
  res = div(res, prod([cont.edge_to_gen[e] for e in cont.edge_vect]))
  res = sum([homogeneous_component(res, i) for i in 0:cont.N_deg-cont.ne])
  res = res * (-1)^count(v -> st.col[v] == 3, 1:nv)
  return res
end

function _pixton_fraction(st::StratumTree)
  cont = st.cont_t[1]
  R = parent(cont.cont_t)
  res = one(R)
  denom = Vector{PolyType}()
  cps = _critical_paths(st)
  nv = Graphs.nv(st.GG)
  for p in cps
    path_sum = sum(vcat([one(R)], [cont.edge_to_gen[e] for e in p]))
    push!(denom, path_sum)
  end
  for v in 1:nv
    st.gen[v] != 0 && continue
    for w in 1:nv
      st.col[w] != 2 && continue # v has genus zero, w is a generalized root.
      p = Graphs.a_star(st.GG, v, w)
      any(e -> st.col[Graphs.src(e)] != 1, p) && continue # Check the path has no interior vertices of positive genus.
      res *= (sum(vcat([one(R)], [cont.edge_to_gen[e] for e in p])))^(Graphs.degree(st.GG, v)-2)
      # res = sum([homogeneous_component(res, i) for i in 0:cont.N_deg-cont.ne])
    end
  end
  res = res // prod(denom)
  res = res * (-1)^count(v -> st.col[v] == 3, 1:nv)
  return res
end

function _pixton_lycka(st::StratumTree)
  nv = Graphs.nv(st.GG)
  # convention: a path goes from lower vertex to higher vertex index.

  # Step 1: generate all the paths between different colored positive genus vertices.
  paths = Vector{Vector{GraphsEdgeType}}()
  for v in 1:nv
    st.col[v] == 1 && continue
    for w in v+1:nv
      st.col[w] == 1 && continue
      st.col[v] == st.col[w] && continue

      p = Graphs.a_star(st.GG, v, w)
      # all interior vertices must have genus zero.
      if length(p) > 1
        any(e -> st.col[Graphs.dst(e)] != 1, p[1:length(p)-1]) && continue
      end
      push!(paths, p)
      println("Added path $p")
    end
  end

  # Step 2: compute all pairwise intersections and count for each path arising that way, what its exponent is.
  pathExps = Dict{Vector{GraphsEdgeType}, Int64}()
  # path_edges = Dict{Vector{Int64}, Vector{GraphsEdgeType}}()
  for path_collection in Combinatorics.combinations(paths)
    k = length(path_collection)
    ip = _intersect_paths(path_collection)
    length(ip) == 0 && continue
    if !haskey(pathExps, ip)
      pathExps[ip] = 0
    end
    pathExps[ip] += (-1)^k
  end
  println(pathExps)

  # Step 3: compute resulting polynomial.
  cont = st.cont_t[1]
  R = parent(cont.cont_t)
  res = one(R)
  denoms = Vector{PolyType}()
  for (p, exp_of_p) in pathExps
    path_factor = sum(vcat([one(R)], [cont.edge_to_gen[e] for e in p]))
    if exp_of_p > 0
      res *= path_factor ^ exp_of_p
    elseif exp_of_p < 0
      push!(denoms, path_factor ^ (-exp_of_p))
    end
  end
  res *= _denominator_to_power_series(denoms, cont.N_deg)
  res = div(res, prod([cont.edge_to_gen[e] for e in cont.edge_vect]))
  res = sum([homogeneous_component(res, i) for i in 0:cont.N_deg-cont.ne])
  #TODO: what is the right sign?
  #res = res * (-1)^count(v -> st.col[v] == 3, 1:nv)
  return res
end

function _intersect_paths(pc::Vector{Vector{GraphsEdgeType}})::Vector{GraphsEdgeType}
  return filter(e -> all(p -> e in p || Graphs.reverse(e) in p, pc), pc[1])
end

