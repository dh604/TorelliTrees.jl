GraphsEdgeType = Graphs.SimpleGraphs.SimpleEdge{Int64}
PolyType = MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}
PolyFracType = AbstractAlgebra.Generic.FracFieldElem{PolyType}

abstract type Abstract_ST_smoothing end

struct Cont_t_container
  cont_t::PolyType
  ne::Int64
  nl::Int64
  N_deg::Int64
  edge_vect::Vector{GraphsEdgeType}
  edge_to_gen::Dict{GraphsEdgeType, PolyType}
end

struct StratumTree
  ls::Vector{Int64}
  GG::Graphs.SimpleGraph{Int64}
  col::Vector{Int64}
  gen::Vector{Int64}
  min_smoothings::Vector{Abstract_ST_smoothing}
  all_smoothings::Vector{Abstract_ST_smoothing}
  cont_t::Vector{Cont_t_container}
end

struct ST_smoothing <: Abstract_ST_smoothing
  ST1::StratumTree # ST1 is a smoothing of ST2
  ST2::StratumTree # ST2 is a degeneration of ST1
  vertex_map::Vector{Vector{Int64}} # vertex_map[v] is the vertices of ST2 corresponding to v in ST1.
  edge_map::Dict{GraphsEdgeType, GraphsEdgeType} # should contain both orientations: E(ST1) -> E(ST2)
end

function Base.show(io::IO, ST::StratumTree)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "Stratum tree")
  else
    # nested printing allowed, preferably terse
    print(io, "Graph $(ST.ls), col $(ST.col), gen $(ST.gen), auts: $(count_auts(ST)), min smoothings: $(length(ST.min_smoothings)), all smoothings: $(length(ST.all_smoothings))")
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", ST::StratumTree)

  print(io, "Graph $(ST.ls), col $(ST.col), gen $(ST.gen), auts: $(count_auts(ST)), min smoothings: $(length(ST.min_smoothings)), all smoothings: $(length(ST.all_smoothings)), codim = $(_cont_deg(ST))")
end

function Base.show(io::IO, sm::ST_smoothing)

  if Oscar.is_terse(io)
    # no nested printing
    print(io, "Stratum tree smoothing")
  else
    # nested printing allowed, preferably terse
    println(io, "Stratum tree smoothing")
    print(io, "  From: ")
    print(io, sm.ST1)
    print(io, "\n    To: ")
    print(io, sm.ST2)
    print(io, "\nV. map: ")
    print(io, sm.vertex_map)
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", sm::ST_smoothing)
  println(io, "Stratum tree smoothing")
  print(io, "      From: ")
  print(io, sm.ST1)
  print(io, "\n        To: ")
  print(io, sm.ST2)
  print(io, "\nVertex map: ")
  print(io, sm.vertex_map)
  print(io, "\n  Edge map: ")
  for (e1, e2) in sm.edge_map
    @req Graphs.reverse(sm.edge_map[e1]) == Graphs.reverse(e2) "Edge map $e1 => $e2 incompatible under Graphs.reverse(...)"
    Graphs.src(e1) > Graphs.src(e2) && continue
    print(io, "$((Graphs.src(e1), Graphs.dst(e1))) => $((Graphs.src(e2), Graphs.dst(e2)))")
    print(io, ", ")
  end
  print(io, ";")
end