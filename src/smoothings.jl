
# Return the set of minimal smoothings of the given stratum tree.
# This assumes that the tree is indeed a valid stratum tree.
function _get_minimal_smoothings(ST::StratumTree)
  minimal_edge_sets = Vector{Vector{GraphsEdgeType}}()
  GG = ST.GG
  col = ST.col
  ls = ST.ls
  for e in Graphs.edges(GG)
    v = Graphs.src(e)
    w = Graphs.dst(e)
    col_v = ST.col[v]
    col_w = ST.col[w]

    # Cannot contract edge connecting two positive genus vertices of different colors.
    col_v != 1 && col_w != 1 && col_v != col_w && continue
    # Contracting an edge between two genus zero vertices is always possible.
    if col_v == 1 && col_w == 1
      push!(minimal_edge_sets, [e])
      continue
    end
    # By assumption that ST is extremal, one color is 1, the other one isn't.
    # Wlog assume col_v == 1 and col_w != 1.
    if col_w == 1
      col_w = col_v
      col_v = 1
      tmp = v
      v = w
      w = tmp
    end
    # From now on, col_v = 1, col_w != 1.
    contracted_edges_for_e = [e]
    col_tmp = copy(col)
    col_tmp[v] = col_w # Possibly split genus zero component of v.
    cpts = _genus_zero_components(ls, col_tmp)

    # Check which components we must contract - in two steps:
    # Step 1: contract all genus zero components with at most one color that we have created by contracting e.
    for (cpt_verts, col_count) in cpts
      col_count >= 2 && continue
      # Contract this component, as it has only one adjacent color.
      for e2 in Graphs.edges(GG)
        (Graphs.src(e2) in cpt_verts || Graphs.dst(e2) in cpt_verts) && push!(contracted_edges_for_e, e2)
      end
    end
    # Step 2: contract all edges incident to v whose color is col_w.
    for u in Graphs.all_neighbors(GG, v)
      u == w && continue
      col[u] == col_w && push!(contracted_edges_for_e, Graphs.Edge(v, u))
    end
    push!(minimal_edge_sets, contracted_edges_for_e)
  end
  # Finally, turn edge lists into vertex partitions.
  function edge_set_to_vertex_partition(edge_set::Vector{GraphsEdgeType})
    return sort(unique(vcat(Graphs.src.(edge_set), Graphs.dst.(edge_set))))
  end
  vertex_sets = [edge_set_to_vertex_partition(edge_set) for edge_set in minimal_edge_sets]
  return sort(unique(vertex_sets))
end

# Return a tuple (Graphs.jl-Graph, col, gen) obtained by smoothing the StratumTree ST at the vertex set sm.
# This assumes that the vertex set determines a subtree.
function _smooth(ST::StratumTree, sm::Vector{Int64})
  nv = length(ST.col)
  nv_sm = nv - length(sm) + 1 # the new vertex corresponding to the contracted vertices has the maximum index.
  
  edge_map = Dict{GraphsEdgeType, GraphsEdgeType}()

  # Track where the vertices go under the contraction map from the degeneration to the smoothing.
  vert_contraction_map = zeros(Int64, nv)
  tmp_ctr = 1
  for v in 1:nv
    if v in sm
      vert_contraction_map[v] = nv_sm
    else
      vert_contraction_map[v] = tmp_ctr
      tmp_ctr += 1
    end
  end
  # Build the contracted tree
  GG_sm = Graphs.SimpleGraph(nv_sm)
  col_sm = zeros(Int64, nv_sm)
  gen_sm = zeros(Int64, nv_sm)
  
  # Step 1: Create edges in smoothened graph
  for e in Graphs.edges(ST.GG)
    v = vert_contraction_map[Graphs.src(e)]
    w = vert_contraction_map[Graphs.dst(e)]
    # discard contracted edges
    v == w && v == nv_sm && continue
    # sanity check:
    @req v != w || v == nv_sum "Requesting a smoothing that causes self-edges at vertices!"
    # create non-contracted edges:
    e_new = Graphs.Edge(v, w)
    Graphs.add_edge!(GG_sm, e_new)
    edge_map[e_new] = e
    edge_map[Graphs.reverse(e_new)] = Graphs.reverse(e)
  end

  # Step 2: Color and genus for contracted vertex:
  col_sm[nv_sm] = 1
  gen_sm[nv_sm] = 0
  for v in sm
    if col_sm[nv_sm] == 1 && ST.col[v] != 1
      col_sm[nv_sm] = ST.col[v]
    end
    gen_sm[nv_sm] += ST.gen[v]
  end
  # Step 3: Color and genus for non-contracted vertices
  for v in 1:nv
    v_new = vert_contraction_map[v]
    v_new == nv_sm && continue
    col_sm[v_new] = ST.col[v]
    gen_sm[v_new] = ST.gen[v]
  end

  # Create vertex map, i.e. for each vertex of GG_sm the set of vertices it corresponds to in GG.
  vertex_map = Vector{Vector{Int64}}()
  for v in 1:nv
    !(v in sm) && push!(vertex_map, [v])
  end
  push!(vertex_map, sm)

  return GG_sm, col_sm, gen_sm, vertex_map, edge_map
end

# Concatenate two smoothings sm1 and sm2, where sm1 degenerates T0 to T1 and sm2 degenerates T1 to T2.
function _concatenate(sm1::ST_smoothing, sm2::ST_smoothing)::ST_smoothing
  @req sm1.ST2 == sm2.ST1 "Canont concatenate incompatible smoothings"
  # Concatenate edge map
  edge_map = Dict{GraphsEdgeType, GraphsEdgeType}()
  for e in keys(sm1.edge_map)
    edge_map[e] = sm2.edge_map[sm1.edge_map[e]]
  end
  # Concatenate vertex map
  vertex_map = Vector{Vector{Int64}}()
  for v0 in 1:Graphs.nv(sm1.ST1.GG)
    vert_v0_part = Int64[]
    for v1 in sm1.vertex_map[v0]
      append!(vert_v0_part, sm2.vertex_map[v1])
    end
    push!(vertex_map, sort(unique(vert_v0_part)))
  end
  return ST_smoothing(sm1.ST1, sm2.ST2, vertex_map, edge_map)
end

# Use this on the result of stratum_trees(...) to calculate all smoothings.
function _add_all_smoothings_from_min_smoothings(STs::Vector{StratumTree})
  for st in STs
    for sm1 in st.min_smoothings
      !_contains_smoothing_up_to_ST1_aut(st.all_smoothings, sm1) && push!(st.all_smoothings, sm1)
      for sm0 in sm1.ST1.all_smoothings
        sm01 = _concatenate(sm0, sm1)
        !_contains_smoothing_up_to_ST1_aut(st.all_smoothings, sm01) && push!(st.all_smoothings, sm01)
      end
    end
  end
end

function _contains_smoothing_up_to_ST1_aut(smoothings::Vector{Abstract_ST_smoothing}, sm::Abstract_ST_smoothing)::Bool
  for sm_prev in smoothings
    if sm_prev.ST1 == sm.ST1 && sm_prev.ST2 == sm.ST2
      ST1 = sm.ST1
      ST1_aut_reln(u, v) = ST1.col[u] == ST1.col[v] && ST1.gen[u] == ST1.gen[v] && sm_prev.vertex_map[u] == sm.vertex_map[v]
      if Graphs.Experimental.has_isomorph(ST1.GG, ST1.GG; vertex_relation = ST1_aut_reln)
        return true
      end
    end
  end
  return false
end

function _count_auts(sm::Abstract_ST_smoothing)::Int64
  
end