# This file is an adaption of adm_interface.jl optimized specifically for the case {2, g-2}.
# Instead of computing actual cycles and products of lambda classes, it restricts itself to
# admcycles' hodge_integral function.

function adm_code_hodge(STs::Vector{StratumTree}, prefix::String, filename::String)
  g_tot = sum(STs[1].gen)
  io = open(filename * ".sage", "w")
  write(io, "from admcycles import *\n\n")
  write(io, "$(prefix)_k1_integral = 0;\n")
  i = 0
  for st in STs
    i += 1
    lns = _contrib_to_adm_hodge(st, prefix, i, length(STs))
    for l in lns
      write(io, l)
      write(io, "\n")
    end
    write(io, "\n")
  end
  write(io, "\nprint('Result for $prefix: ')\n")
  write(io, "print($(prefix)_k1_integral);\n\n")
  write(io, "# Save result in text file:\n")
  write(io, "f = open('$(prefix)_k1_integral.txt', 'w');\n")
  write(io, "f.write(str($(prefix)_k1_integral));\n")
  write(io, "f.close();\n")
  close(io)
end

function _contrib_to_adm_hodge(st::StratumTree, prefix::String, treeNr::Int64, nTrees::Int64; print_compact::Bool=false, multithreading::Bool=false, thread_nr::Int64=-1)
  cont = st.cont_t[1]
  cont_t = cont.cont_t
  GG = st.GG
  nv = Graphs.nv(GG)
  ne = Graphs.ne(GG)
  g_tot = sum(st.gen)
  n_vars = g_tot + 2*ne

  # Prepare variable names for debugging
  variable_names = Vector{String}(undef, n_vars)
  adm_names = Vector{String}(undef, n_vars)

  # Assign edges to half edges and half edges to vertices
  edge_to_half_edge = Dict{GraphsEdgeType, Int64}() # Assign the index of the corresponding half edge to an oriented edge
  half_edge_to_vertex = Vector{Int64}() # Assign the vertex to a half-edge.
  half_edge_index = 1
  for v in 1:nv
    half_edges_at_v = Vector{Int64}()
    for w in Graphs.all_neighbors(GG, v)
      push!(half_edges_at_v, half_edge_index)
      edge_to_half_edge[Graphs.Edge(v, w)] = half_edge_index
      push!(half_edge_to_vertex, v)
      half_edge_index += 1
    end
  end
  vertex_count = zeros(Int64, nv)
  psi_gen_to_vertex_count = zeros(Int64, 2*ne)
  for i in 1:2*ne
    v = half_edge_to_vertex[i]
    vertex_count[v] += 1
    g = st.gen[v]
    n = Graphs.degree(GG, v)
    adm_names[i] = "psiclass($(vertex_count[v]), $g, $n)"
    variable_names[i] = "ψ_{$v,$i}"
    psi_gen_to_vertex_count[i] = vertex_count[v]
  end

  # Assign each generator of T to a vertex of GG:
  # First the half edges in order of half_edge_to_vertex, then the lambda classes
  gen_to_vertex = zeros(Int64, n_vars)
  gen_degrees = ones(Int64, n_vars)
  gen_to_vertex[1:2*ne] = half_edge_to_vertex
  i = 1
  for v in 1:nv
    g = st.gen[v]
    g == 0 && continue
    n = Graphs.degree(GG, v)
    gen_to_vertex[2*ne+i:2*ne+i+g-1] = [v for _ in 1:g]
    gen_degrees[2*ne+i:2*ne+i+g-1] = 1:g
    variable_names[2*ne+i:2*ne+i+g-1] = ["λ_{$v,$i}" for i in 1:g]
    adm_names[2*ne+i:2*ne+i+g-1] = ["lambdaclass($i, $g, $n)" for i in 1:g]
    i += g
  end

  # Prepare substitution: First edges, then l's, then ci's.
  T, gens = graded_polynomial_ring(QQ, variable_names, gen_degrees)
  ADM_T, ADM_gens = graded_polynomial_ring(QQ, adm_names, gen_degrees)
  edge_subst::Vector{PolyType} = [-gens[edge_to_half_edge[e]] - gens[edge_to_half_edge[Graphs.reverse(e)]] for e in cont.edge_vect]
  l_subst::Vector{PolyType} = zeros(T, cont.nl)
  ccl = _total_chern_class(st, nv, gens, gen_to_vertex)
  c_subst::Vector{PolyType} = [homogeneous_component(ccl, i) for i in 1:cont.N_deg]
  subst::Vector{PolyType} = vcat(edge_subst, l_subst, c_subst)
  F = evaluate(cont_t, subst)
  # println(F)

  # Filter out terms that do not contribute for degree reasons.
  F_filtered = zero(F)
  for t in terms(F)
    e = collect(exponents(t))[1]
    filtered_out = false
    for v in 1:nv
      e_v = [gen_to_vertex[i] == v ? e[i] : 0 for i in 1:n_vars]
      d_v = sum(e_v .* gen_degrees)
      if d_v > 2*st.gen[v] - 3 + Graphs.degree(GG, v)
        # println("Filtered out $t for degree $d_v at vertex $v.")
        filtered_out = true
        break
      end
    end
    if !filtered_out
      F_filtered += t
    end
  end
  print_compact && println(F_filtered)

  # Finally, convert result to admcycles code.
  output_lines = Vector{String}()
  push!(output_lines, "$(prefix)_T$(treeNr)_term = 0;")
  push!(output_lines, "")
  lt = length(terms(F_filtered))
  tctr = 0
  for t in terms(F_filtered)
    tctr += 1
    push!(output_lines, string("$(prefix)_T$(treeNr)_term += $(_adm_hodge_integral_product(st, t, gen_to_vertex, ADM_gens, gen_degrees, psi_gen_to_vertex_count, variable_names));"))
    push!(output_lines, "print('Computed term $tctr / $lt of tree $treeNr / $nTrees');")
  end
  push!(output_lines, "$(prefix)_k1_integral += 1/$(count_auts(st)) * $(prefix)_T$(treeNr)_term;")
  return output_lines
end

function _adm_hodge_integral_product(st::StratumTree, t::PolyType, gen_to_vertex::Vector{Int64}, ADM_gens, gen_degrees, psi_gen_to_vertex_count, variable_names)::String
  n_vars = length(ADM_gens)
  GG = st.GG
  nv = Graphs.nv(GG)
  e = collect(exponents(t))[1]
  res = ""
  for v in 1:nv
    g = st.gen[v]
    n = Graphs.degree(GG, v)
    e_v = [gen_to_vertex[i] == v ? e[i] : 0 for i in 1:n_vars]
    d_v = sum(e_v .* gen_degrees)
    if v == 1
      res *= string(coeff(t, monomial(t, 1)))
      res *= " * "
    end

    # Turn monomial at v into hodge integral.

    lambdas = [g] #  we use the lambda_g evaluation.
    psis = zeros(Int64, Graphs.degree(GG, v))
    for i in 1:n_vars

      gen_to_vertex[i] != v && continue

      if variable_names[i][1] == 'λ'
        append!(lambdas, gen_degrees[i] .* ones(Int64, e_v[i]))
      elseif variable_names[i][1] == 'ψ'
        vert_count = psi_gen_to_vertex_count[i]
        psis[vert_count] += e_v[i]
      end
    end
    
    if d_v == 2*st.gen[v] - 4 + Graphs.degree(GG, v)
      kappas = [1]
    else
      kappas = Int64[]
    end

    res *= "hodge_integral($g, $(_int_string_to_admcycles(lambdas)), $(_int_string_to_admcycles(psis)), $(_int_string_to_admcycles(kappas)), n=$n)"
    if v != nv
      res *= " * "
    end
  end
  return res * ""
end

function _int_string_to_admcycles(s::Vector{Int64})::String
  if length(s) == 0
    return "[]"
  else
    return "$s"
  end
end