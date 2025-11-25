
# Brute force version for only taking one representative per isomorphism class of coloured tree.
function stratum_trees(g::Vector{Int64}; printStuff::Bool=false, draw::Bool=false, folder_name::String="", compute_smoothings::Bool=true)
  @req all(i -> i > 0, g) "g must have strictly positive entries"
  g_tot = sum(g)
  codim = div(g_tot*(g_tot+1) - sum(gk -> gk*(gk+1), g), 2)
  printStuff && println("Codimension $codim")
  max_ne = codim
  max_nv = max_ne + 1 # since all stable graphs are trees.
  n_cols = length(g)
  min_per_cols = vcat([-1], ones(Int64, n_cols))
  max_per_cols = vcat([-1], g)

  return _stratum_trees(max_nv, g, min_per_cols, max_per_cols; draw=draw, folder_name=folder_name, printStuff=printStuff, compute_smoothings=compute_smoothings)
end

# Brute force version for only taking one representative per isomorphism class of coloured tree.
function irred_stratum_trees(g::Vector{Int64}; draw::Bool=false, folder_name::String="", printStuff::Bool=false, compute_smoothings::Bool=true)
  @req all(i -> i > 0, g) "g must have strictly positive entries"
  g_tot = sum(g)
  codim = div(g_tot*(g_tot+1) - sum(gk -> gk*(gk+1), g), 2)
  printStuff && println("Codimension $codim")
  max_ne = codim
  max_nv = max_ne + 1 # since all stable graphs are trees.
  n_cols = length(g)
  min_per_cols = vcat([-1], ones(Int64, n_cols))
  max_per_cols = vcat([0], g)

  return _stratum_trees(max_nv, g, min_per_cols, max_per_cols; draw=draw, folder_name=folder_name, printStuff=printStuff, compute_smoothings=compute_smoothings)
end

# Brute force version for only taking one representative per isomorphism class of coloured tree.
function _stratum_trees(max_nv::Int64, g::Vector{Int64}, min_per_cols::Vector{Int64}=Int64[], max_per_cols::Vector{Int64}=Int64[]; printStuff::Bool=false, draw::Bool=false, folder_name::String="", compute_smoothings::Bool=true)

  n_cols = length(g)
  @req all(i -> i > 0, g) "g must have strictly positive entries"

  rand_col() = RGB(rand()*0.7+0.1, rand()*0.7+0.1, rand()*0.1+0.1)
  if draw
    fixed_cols = [RGB(0.8,0.8,0.8), RGB(0.85, 0.65, 0.13), RGB(0, 0.5, 0.5), RGB(0.5, 0.8, 0.5)]
    cols = vcat(fixed_cols, [rand_col() for _ in 1:n_cols-length(fixed_cols)+1])
  # println(cols)
  end

  ans::Int64 = 0
  res = Vector{StratumTree}()

  isempty(max_per_cols) && (max_per_cols = [-1 for _ in 1:n_cols+1]) # -1 means no maximum is imposed.
  isempty(min_per_cols) && (min_per_cols = [-1 for _ in 1:n_cols+1]) # -1 means no minimum is imposed.

  prev_st_for_ls = Dict{Vector{Int64}, Vector{StratumTree}}()

  for ls in Iterators.flatten([TreeIt(i) for i in 2:max_nv]) # generate the level sequences
    
    ls = copy(ls) # This is crucial as the iterator subsequently changes the level sequence object.

    G = LStoGraph(ls)
    GG = Graphs.SimpleGraph(Oscar.adjacency_matrix(G))

    ls_cols = Vector{Vector{Int64}}()

    CI = _extremal_colourings_with_repetitions(ls, GG, n_cols, min_per_cols, max_per_cols) # generation of colorings (with repeated isomorphism classes)
    for col in CI
      col = copy(col) # To be on the safe side.
      # Only consider one coloring per isomorphism class.
      any(col_old -> _is_isomorphic_coloring(GG, col_old, col), ls_cols) && continue
      push!(ls_cols, col)

      ls_col_gens = Vector{Vector{Int64}}()
      for gen in _stratum_trees_with_repetitions(GG, col, g)
        gen = copy(gen) # To be on the safe side.

        # Only consider one per isom class.
        any(old_gen -> _is_isomorhic_gen(GG, col, gen, old_gen), ls_col_gens) && continue
        push!(ls_col_gens, gen)

        # Found a new stratum tree!
        st = StratumTree(ls, GG, col, gen, Vector{ST_smoothing}(), Vector{ST_smoothing}(), Vector{Cont_t_container}(undef, 1))
        # println("\nCreating tree $st with edges:")
        # println(collect(Graphs.edges(GG)))
        # println()
        ans += 1
        push!(res, st)
        # Add stratum tree to list, indexed by level sequence.
        if haskey(prev_st_for_ls, ls)
          push!(prev_st_for_ls[ls], st)
        else
          prev_st_for_ls[ls] = [st]
        end
        # println("prev_st_for_ls[$ls] = $(prev_st_for_ls[ls])")

        auts = count_auts(st)
        printStuff && println("\nGraph $ls, col $col, gen $gen and $auts auts.")
        
        if compute_smoothings
          min_smoothings = _get_minimal_smoothings(st)
          printStuff && println("Minimal smoothings: $(min_smoothings)")

          # For each minimum smoothing, create the corresponding object.
          for sm in min_smoothings
            # println("smooth $sm")
            # Step 1: create smoothed graph
            GG_sm, col_sm, gen_sm, vertex_map, edge_map = _smooth(st, sm)
            nv_sm = Graphs.nv(GG_sm)
            vertex_map_after_iso = Vector{Vector{Int64}}(undef, nv_sm)
            edge_map_after_iso = Dict{GraphsEdgeType, GraphsEdgeType}()

            # println("GG_sm:")
            # println("nv_sm = $nv_sm")
            # println(collect(Graphs.edges(GG_sm)))
            # println("gen_sm = $gen_sm")
            # println("col_sm = $col_sm")


            # Step 2: find isomorphism type of tree in list of previous StratumTrees.
            ls_sm = Vector{Int64}()
            # println("reset ls_sm")
            for (ls_candidate, st_candidates) in prev_st_for_ls
              st_candidate = st_candidates[1]
              if Graphs.Experimental.has_isomorph(st_candidate.GG, GG_sm)
                ls_sm = ls_candidate
                # println("Approved tree candidate: $st_candidate, which has edges:")   
                # println(collect(Graphs.edges(st_candidate.GG)))
                # println("ls_candidate = $ls_candidate, st_candidate.ls = $(st_candidate.ls)")
                break
              end
            end
            @req length(ls_sm) > 0 "Smoothing not found in isomorphism types of previous trees!"
            # println("ls_sm = $ls_sm")

            # Step 3: Find full stratum tree corresponding to smoothing.
            st_sm::Union{Nothing, StratumTree} = nothing
            for st_candidate in prev_st_for_ls[ls_sm]
              if _is_isomorphic_col_and_gen(st_candidate.GG, GG_sm, st_candidate.col, col_sm, st_candidate.gen, gen_sm)
                st_sm = st_candidate
                break
              end
            end
            @req !isnothing(st_sm) "Smoothing not found in previous stratum trees!"

            # Step 4: use isomorphism to rewrite vertex_map, edge_map.
            iso = _get_isom(st_sm.GG, GG_sm, st_sm.col, col_sm, st_sm.gen, gen_sm)
            iso_vert_map = zeros(Int64, nv_sm)
            for (a, b) in iso
              iso_vert_map[a] = b
            end
            for v in 1:nv_sm
              vertex_map_after_iso[v] = vertex_map[iso_vert_map[v]]
            end
            for e in Graphs.edges(st_sm.GG)
              v = Graphs.src(e)
              w = Graphs.dst(e)
              e_middle = Graphs.Edge(iso_vert_map[v], iso_vert_map[w])
              e_target = edge_map[e_middle]
              edge_map_after_iso[e] = e_target
              edge_map_after_iso[Graphs.reverse(e)] = Graphs.reverse(e_target)
            end

            final_sm_object = ST_smoothing(st_sm, st, vertex_map_after_iso, edge_map_after_iso)
            push!(st.min_smoothings, final_sm_object)
          end
        end

        draw && draw_pdf(st, string("tree_", ans, ".pdf"), 100, 100, cols, auts, compute_smoothings)
      end
    end
  end
  if draw
    mod5 = 0
    for i in 1:ans
      println(string("\\includegraphics[width=3.2cm]{", folder_name, "tree_", i, ".pdf}"))
      mod5 += 1
      if mod5 == 5
        mod5 = 0
        println()
      end
    end
  end
  if compute_smoothings
    _add_all_smoothings_from_min_smoothings(res)
  end
  return res
end

function count_auts(st::StratumTree)
  col_and_gen_rel(u, v) = (st.col[u] == st.col[v]) && (st.gen[u] == st.gen[v])
  return Graphs.Experimental.count_isomorph(st.GG, st.GG; vertex_relation = col_and_gen_rel)
end

function _get_isom(GG1::Graphs.SimpleGraph{Int64}, GG2::Graphs.SimpleGraph{Int64}, col1::Vector{Int64}, col2::Vector{Int64}, gen1::Vector{Int64}, gen2::Vector{Int64})::Vector{Tuple{Int64, Int64}}
  color_rel(u, v) = col1[u] == col2[v] && gen1[u] == gen2[v]
  isom_iterator = Graphs.Experimental.all_isomorph(GG1, GG2; vertex_relation = color_rel)
  for iso in isom_iterator
    return iso
  end
end

# Check if the two colorings including genus distributions are isomorphic.
function _is_isomorphic_col_and_gen(GG1::Graphs.SimpleGraph{Int64}, GG2::Graphs.SimpleGraph{Int64}, col1::Vector{Int64}, col2::Vector{Int64}, gen1::Vector{Int64}, gen2::Vector{Int64})::Bool
  color_rel(u, v) = col1[u] == col2[v] && gen1[u] == gen2[v]
  return Graphs.Experimental.has_isomorph(GG1, GG2; vertex_relation = color_rel)
end

# Check if the two colorings col1 and col2 on GG are isomorphic.
function _is_isomorphic_coloring(GG::Graphs.SimpleGraph{Int64}, col1::Vector{Int64}, col2::Vector{Int64})::Bool
  color_rel(u, v) = col1[u] == col2[v]
  return Graphs.Experimental.has_isomorph(GG, GG; vertex_relation = color_rel)
end

# Check if the two colorings col1 and col2 on GG are isomorphic.
function _is_isomorhic_gen(GG::Graphs.SimpleGraph{Int64}, col::Vector{Int64}, gen1::Vector{Int64}, gen2::Vector{Int64})::Bool
  color_rel(u, v) = col[u] == col[v] && gen1[u] == gen2[v]
  return Graphs.Experimental.has_isomorph(GG, GG; vertex_relation = color_rel)
end

# Here col_sums has length n_cols, not n_cols+1.
function _stratum_trees_with_repetitions(GG::Graphs.SimpleGraph{Int64}, col::Vector{Int64}, col_sums::Vector{Int64})
  n_cols = length(col_sums)
  nv = length(col)
  cts = [count(i -> i == c, col) for c in 2:n_cols+1] # cts[c] = number of vertices on GG with color c.

  # Iterate over correct compositions
  comp_iterator = Iterators.product((Oscar.compositions(col_sums[i], cts[i]) for i in 1:n_cols)...,)
  
  # println("col=$col")
  # Turn composition into genus list.
  function comp_to_gen(comp)
    res = zeros(Int64, nv)
    col_counter = ones(Int64, n_cols)
    for v in 1:nv
      col_v = col[v]
      if col_v == 1
        res[v] = 0
        continue
      end
      # println("comp: $comp, col_v = $col_v, col_counter[col_v-1]=$(col_counter[col_v-1]), res=$res")
      res[v] = comp[col_v-1][col_counter[col_v-1]]
      col_counter[col_v-1] += 1
    end
    return res
  end

  return Iterators.map(comp_to_gen, comp_iterator)
end

# Return iterator over all colorings which are extremal and permissible.
# This may return multiple representatives of the same isomorphism type of colored tree.
function _extremal_colourings_with_repetitions(ls::Vector{Int64}, GG::Graphs.SimpleGraph{Int64}, n_cols::Int64, min_per_cols::Vector{Int64}, max_per_cols::Vector{Int64})
  nv = Graphs.nv(GG)
  col_distributions = _col_distribution_iterator(nv, n_cols, min_per_cols, max_per_cols)
  return Iterators.filter(col -> _is_permissible_coloring(GG, col) && _is_extremal(ls, col), col_distributions)
end

# Return iterator over all colorings within the specified ranges, regardless of whether it is extremal and permissible.
function _col_distribution_iterator(nv::Int64, n_cols::Int64, min_per_cols::Vector{Int64}, max_per_cols::Vector{Int64})
  colour_compositions = _permissible_colour_compositions(nv, n_cols, min_per_cols, max_per_cols) # iterator over permissible [col1, col2, ...] such that col1+col2+...=nv
  composition_to_vector = cts -> vcat(([i for _ in 1:cts[i]] for i in 1:n_cols+1)...,) # [col1 times 1, col2 times 2, ...]
  composition_to_iterator = cts -> Combinatorics.multiset_permutations(composition_to_vector(cts), nv) # Iterating over permutations thereof
  col_distributions = Iterators.flatten(map(cts -> composition_to_iterator(cts), colour_compositions))
  return col_distributions
end


# return iterator over [col1, col2, ...] within the specified range such that col1+col2+...=nv
function _permissible_colour_compositions(nv::Int64, n_cols::Int64, min_per_cols::Vector{Int64}, max_per_cols::Vector{Int64})
  all_compositions = Oscar.weak_compositions(nv, n_cols+1)
  max_check = cts -> all(i -> max_per_cols[i] < 0 || max_per_cols[i] >= cts[i], 1:n_cols+1)
  min_check = cts -> all(i -> min_per_cols[i] < 0 || min_per_cols[i] <= cts[i], 1:n_cols+1)
  Iterators.filter(cts -> max_check(cts) && min_check(cts), all_compositions)
end

# Check if no two vertices of colour != 1 are next to each other.
function _is_permissible_coloring(GG::Graphs.SimpleGraph{Int64}, col::Vector{Int64})::Bool
  for e in Graphs.edges(GG)
    v = Graphs.src(e)
    w = Graphs.dst(e)
    col[v] != 1 && col[w] != 1 && col[v] == col[w] && return false
  end
  return true
end

# Check if every color 1 component is adjacent to at least two different colours and every color 1 vertex has degree at least three.
function _is_extremal(ls::Vector{Int64}, col::Union{Vector{Int64}, Tuple{Vararg{Int64}}})::Bool
  nv = length(ls)
  G = LStoGraph(ls)
  for v in 1:nv
    col[v] == 1 && degree(G, v) < 3 && return false
  end
  g0cpts = _genus_zero_components(ls, col; G=G)
  return all(cpt -> cpt[2] >= 2, g0cpts)
end

# Return a vector whose entries are the genus zero components (given as vertex set) and the number of adjacent colors.
function _genus_zero_components(ls::Vector{Int64}, col::Union{Vector{Int64}, Tuple{Vararg{Int64}}}; G::Union{Oscar.Graph, Nothing}=nothing)
  res = Vector{Tuple{Vector{Int64}, Int64}}()
  nv = length(ls)
  if isnothing(G)
    G = LStoGraph(ls)
  end
  adj_cols = [Int64[] for _ in 1:nv]
  to_be_removed = Vector{Edge}()
  for v in 1:nv
    col[v] == 1 && continue
    for w in all_neighbors(G, v)
      col[w] != 1 && continue
      # Now col[v] != 1 and col[w] == 1 holds.
      !(col[v] in adj_cols[w]) && append!(adj_cols[w], col[v])
      push!(to_be_removed, Edge(v, w))
    end
  end
  foreach(e -> Oscar.rem_edge!(G, e), to_be_removed)
  # println("Cut graph: $G")
  # println("Adjacent colors: $adj_cols")
  # check connected components.
  cpts = Oscar.connected_components(G)
  for cpt in cpts
    col[cpt[1]] != 1 && continue # only consider components of color 1.
    adj_cols_in_cpt = unique(Iterators.flatten([adj_cols[w] for w in cpt])) # get all colors adjacent to the component
    push!(res, (cpt, length(adj_cols_in_cpt)))
  end
  return res
end