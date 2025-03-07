export generalized_flag
export generalized_flag2
export generalized_gkm_flag


@doc raw"""
    generalized_gkm_flag(R::RootSystem; S::Vector{RootSpaceElem}=RootSpaceElem[]) -> AbstractGKM_graph{QQFieldElem}

Given a root system ``R`` and a subset ``S`` of the set of simple roots, it construct  the 
generalized flag variety ``G/P``. Here ``G`` is the simple connected complex Lie group 
with root system ``R``, and ``P`` is the parabolic subgroup with root system ``S``.

!!! warning
    Computing this function with root systems of type ``E6``, ``E7``, and ``E8`` may be very slow.

# Examples
```jldoctest
julia> R = root_system([(:A, 1), (:G, 2)])
Root system of rank 3
  of type A1 x G2

julia> S = simple_roots(R)[3];

julia> gp = generalized_gkm_flag(R, [S]);

julia> rank_torus(gp)
5

```
"""
function generalized_gkm_flag(R::RootSystem, S::Vector{RootSpaceElem}=RootSpaceElem[])

  (fams, ordering) = root_system_type_with_ordering(R)

  #TODO add test for different ordering

  W = weyl_group(R) # Weyl group of the root system
  WP = [one(W)] # embedding of W_P into W, called WP

  if !isempty(S)

    #define the Weyl group W_P of the subroot system
    s = simple_roots(R);
    indeces_of_S = findall(j -> s[j] in S, eachindex(s))
    cartan_submatrix = cartan_matrix(R)[indeces_of_S, indeces_of_S]

    #embedding of W_P into W, called WP
    WP = [prod(i -> reflection(s[indeces_of_S[i]]), word(a); init = one(W)) for a in weyl_group(cartan_submatrix)]
  end

  ## Construct the cosets
  # cosets = unique(Set.([b .* WP for b in W]))

  cosets = [WP for _ in 1:div(order(W), length(WP))]
  reprs = [one(W) for _ in 1:length(cosets)]
  index = 1
  for b in W
    any(i -> b in cosets[i], 1:index) && continue
    index += 1
    cosets[index] = b .* cosets[index]
    reprs[index] = reduce((x, y) -> length(x) <= length(y) ? x : y,  cosets[index])
    index == length(cosets) && break
  end

  labs = [replace(repr(r), " " => "") for r in reprs]# repr.(reprs)
  g = Graph{Undirected}(length(labs))
  M = free_module(QQ, mapreduce(_dimension_ambient, +, fams))
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{QQFieldElem}}()

  gen_matrix = AbstractAlgebra.perm(ordering)*block_diagonal_matrix([_generator_matrix(fam) for fam in fams])

  for i in 1:length(reprs)
    r1 = reprs[i]
    for t in positive_roots(R)
      new_rep = reflection(t)*r1
      j = findfirst(index -> (index > i) && (new_rep in cosets[index]), 1:length(reprs)) #change index != i if you want double ways 
      j === nothing && continue
      add_edge!(g, j, i)
      vec = coefficients(t)*gen_matrix
      W[Edge(j,i)] = (-1)*sum(i -> vec[i]*gens(M)[i], 1:rank(M))
    end
  end

  return gkm_graph(g, labs, M, W)
  # return W

end

function generalized_gkm_flag(R::RootSystem, S::Vector{Int64}=Int64[])
end

function generalized_flag2(rt::RootSystem, S::Vector{RootSpaceElem}=RootSpaceElem[])

  (fams, ordering) = root_system_type_with_ordering(rt)

  #TODO add test for different ordering

  W = weyl_group(rt) # Weyl group of the root system
  WP = [one(W)] # embedding of W_P into W, called WP

  if !isempty(S)

    #define the Weyl group W_P of the subroot system
    s = simple_roots(rt);
    indeces_of_S = findall(j -> s[j] in S, eachindex(s))
    cartan_submatrix = cartan_matrix(rt)[indeces_of_S, indeces_of_S]

    #embedding of W_P into W, called WP
    WP = [prod(i -> reflection(s[indeces_of_S[i]]), word(a); init = one(W)) for a in weyl_group(cartan_submatrix)]
  end

  ## Construct the cosets
  # cosets = unique(Set.([b .* WP for b in W]))

  # cosets = [WP for _ in 1:div(order(W), length(WP))]
  reprs = [one(W) for _ in 1:div(order(W), length(WP))]
  index = 1
  for b in W
    any(i -> b^(-1)*reprs[i] in WP, 1:index) && continue
    index += 1
    # cosets[index] = b .* cosets[index]
    
    reprs[index] = b*reduce((x, y) -> length(b*x) <= length(b*y) ? x : y,  WP)
    index == length(reprs) && break
  end

  labs = [replace(repr(r), " " => "") for r in reprs]# repr.(reprs)
  g = Graph{Undirected}(length(labs))
  M = free_module(QQ, mapreduce(_dimension_ambient, +, fams))
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{QQFieldElem}}()

  gen_matrix = block_diagonal_matrix([_generator_matrix(fam) for fam in fams])
  if any(i -> ordering[i] != i, eachindex(ordering))
    gen_matrix = sub(gen_matrix, ordering, 1:number_of_columns(gen_matrix))
  end
#gen_matrix = perm(ordering)*block_diagonal_matrix([_generator_matrix(fam) for fam in fams])


  for i in 1:length(reprs)
    r1 = reprs[i]
    for t in positive_roots(rt)
      # invol = reflection(t)
      j = findfirst(index -> ((index > i) && ((reprs[index])^(-1))*reflection(t)*r1 in WP), 1:length(reprs)) #change index != i if you want double ways 
      j === nothing && continue
      add_edge!(g, j, i)
      vec = coefficients(t)*gen_matrix
      # vec = sub(sub(coefficients(t), [1], ordering)*gen_matrix, [1], Vector(perm(ordering)^(-1), rank(M)))
      W[Edge(j,i)] = (-1)*sum(i -> vec[i]*gens(M)[i], 1:rank(M))
    end
  end
  # println(labs)
  # return (g, labs) #D
  return W

end

function _dimension_ambient(RT::Tuple{Symbol, Int64})::Int64
  # RT = root_system_type(R)[1]

  if RT[1] in (:A, :G)
    return RT[2] + 1
  elseif RT[1] in (:B, :C, :D, :F)
    return RT[2]
  end

  return 8 # RT[1] == :E

end

function _generator_matrix(RT::Tuple{Symbol, Int64})::QQMatrix

  # RT = root_system_type(R)[1]
  n_rows = RT[2]
  n_cols = _dimension_ambient(RT)

  if RT[1] == :E # following Hum75 convention
    M = zero_matrix(QQ, n_rows, n_cols)
    foreach(i -> M[1, i] = (i==1 || i==8) ? QQ(1)//QQ(2) : QQ(-1)//QQ(2), 1:n_cols)
    M[2, 1] = QQ(1)
    M[2, 2] = QQ(1)
    for i in 3:n_rows
      M[i, i-2] = QQ(-1)
      M[i, i-1] = QQ(1)
    end
    return M
  end

  M = diagonal_matrix(QQ(1), n_rows, n_cols)

  if RT[1] == :G  # following Hum75 convention
    M[1, 2] = QQ(-1)
    M[2, 1] = QQ(-2)
    M[2, 3] = QQ(1)
    return M
  end

  for i in 2:n_cols
    M[i-1, i] = QQ(-1)
  end

  if RT[1] in (:A, :B) # # following Hum75 & Wikipedia convention
  else
    if RT[1] == :C # following Hum75 & Wikipedia convention
      M[n_rows, n_cols] = QQ(2)
    elseif RT[1] == :D # following Hum75 & Wikipedia convention
      M[n_rows, n_cols-1] = QQ(1)
    elseif RT[1] == :F # following Wikipedia convention
      M[3, 4] = QQ(0)
      foreach(i -> M[4, i] = QQ(-1)//QQ(2), 1:4)
    end
  end

  return M
end