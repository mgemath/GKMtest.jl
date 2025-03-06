export generalized_flag
export generalized_flag2
# export _generator_matrix

function generalized_flag(rt::RootSystem, S::Vector{RootSpaceElem}=RootSpaceElem[])

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

  ## Find a representative of minimal length for each coset 
  # reprs = [reduce((x, y) -> length(x) <= length(y) ? x : y,  c) for c in cosets]

  ## Involutions of W
  # involutions = [reflection(t) for t in positive_roots(rt)] # alternatives [b for b in W if (order(b) == 2)], [b for b in W if (order(b) == 2 && !(b in WP))]

  # D = Dict{WeylGroupElem, Vector{Tuple{WeylGroupElem, WeylGroupElem}}}()
  # for i in 1:length(reprs)
  #   r1 = reprs[i]
  #   D[r1] = Tuple{WeylGroupElem, WeylGroupElem}[]
  #   for invol in involutions
  #     j = findfirst(index -> (index > i) && (r1*invol in cosets[index]), 1:length(reprs)) #change index != i if you want double ways 
  #     j === nothing && continue
  #     r2 = reprs[j]
  #     push!(D[r1], (r2, invol))
  #   end
  # end
  # println(D)

  labs = [replace(repr(r), " " => "") for r in reprs]# repr.(reprs)
  g = Graph{Undirected}(length(labs))
  M = free_module(QQ, mapreduce(_dimension_ambient, +, fams))
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{QQFieldElem}}()

  gen_matrix = block_diagonal_matrix([_generator_matrix(fam) for fam in fams])
  if any(i -> ordering[i] != i, eachindex(ordering))
    gen_matrix = sub(gen_matrix, ordering, 1:number_of_columns(gen_matrix))
  end
  # gen_matrix = sub(block_diagonal_matrix([_generator_matrix(fam) for fam in fams]), ordering, ordering)

  for i in 1:length(reprs)
    r1 = reprs[i]
    # D[r1] = Tuple{WeylGroupElem, WeylGroupElem}[]
    for t in positive_roots(rt)
      # invol = reflection(t)
      j = findfirst(index -> (index > i) && (reflection(t)*r1 in cosets[index]), 1:length(reprs)) #change index != i if you want double ways 
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

  ## Find a representative of minimal length for each coset 
  # reprs = [reduce((x, y) -> length(x) <= length(y) ? x : y,  c) for c in cosets]

  ## Involutions of W
  # involutions = [reflection(t) for t in positive_roots(rt)] # alternatives [b for b in W if (order(b) == 2)], [b for b in W if (order(b) == 2 && !(b in WP))]

  # D = Dict{WeylGroupElem, Vector{Tuple{WeylGroupElem, WeylGroupElem}}}()
  # for i in 1:length(reprs)
  #   r1 = reprs[i]
  #   D[r1] = Tuple{WeylGroupElem, WeylGroupElem}[]
  #   for invol in involutions
  #     j = findfirst(index -> (index > i) && (r1*invol in cosets[index]), 1:length(reprs)) #change index != i if you want double ways 
  #     j === nothing && continue
  #     r2 = reprs[j]
  #     push!(D[r1], (r2, invol))
  #   end
  # end
  # println(D)

  labs = [replace(repr(r), " " => "") for r in reprs]# repr.(reprs)
  g = Graph{Undirected}(length(labs))
  M = free_module(QQ, mapreduce(_dimension_ambient, +, fams))
  W = Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{QQFieldElem}}()

  gen_matrix = block_diagonal_matrix([_generator_matrix(fam) for fam in fams])
  if any(i -> ordering[i] != i, eachindex(ordering))
    gen_matrix = sub(gen_matrix, ordering, 1:number_of_columns(gen_matrix))
  end

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