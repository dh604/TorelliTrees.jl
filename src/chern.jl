
# Let E and F be total Chern classes of vector bundles (grading comes from the polynomial ring grading.)
# This returns the total chern class of the tensor product of E and F.
function _c_tot_of_product(E::PolyType, F::PolyType)::PolyType
  R = parent(E)
  r = _deg(E)
  q = _deg(F)
  Lambda = zero_matrix(R, q, q)
  Lambda[:,1] = [homogeneous_component(F, i) for i in 1:q]
  foreach(i -> Lambda[i,i+1] = -one(R), 1:q-1)
  # println(Lambda)
  I = identity_matrix(R, q)
  f = I + Lambda
  tmp = identity_matrix(R, q)
  M = zero_matrix(R, q, q)
  for k in 0:r
    M += homogeneous_component(E, r-k) .* tmp
    tmp *= f
  end
  return det(M)
end

function _deg(f::PolyType)::Int64
  @req !iszero(f) "f must be nonzero"
  deg = -1
  for t in terms(f)
    deg = max(degree(t)[1], deg)
  end
  return deg
end

function _total_chern_class(st::StratumTree, nv::Int64, gens::Vector{PolyType}, gen_to_vertex::Vector{Int64})::PolyType
  T = parent(gens[1])
  res = one(T)
  for v in 1:nv
    st.col[v] == 1 && continue
    for w in v+1:nv
      st.col[w] == 1 && continue
      st.col[v] == st.col[w] && continue
      g_v = st.gen[v]
      v_last = findlast(u -> u == v, gen_to_vertex)
      v_indices = v_last-g_v+1:v_last
      g_w = st.gen[w]
      w_last = findlast(u -> u == w, gen_to_vertex)
      w_indices = w_last-g_w+1:w_last
      c_tot_v = one(T) + sum([(-1)^i for i in 1:g_v] .* gens[v_indices])
      c_tot_w = one(T) + sum([(-1)^i for i in 1:g_w] .* gens[w_indices])
      res *= _c_tot_of_product(c_tot_v, c_tot_w)
    end
  end
  return res
end