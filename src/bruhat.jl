export get_bruhat_order_of_generalized_flag, generalized_gkm_schubert

struct BruhatOrder
  order::Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}}
  reversed::Bool # true for descending order
  labels::Vector{String}
end

function _add_to_order(G::AbstractGKM_graph, order::Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}}, elem::Int64, increment::Int64)

    ## Create subgraph for elem
    elem_in_dict = (elem, count('s', G.labels[elem]))
    order[elem_in_dict] = Tuple{Int64, Int64}[]
    
    
    for v in neighbors(G.g, elem)
      if count('s', G.labels[v]) == count('s', G.labels[elem]) + increment 
        push!(order[elem_in_dict], (v, count('s', G.labels[v])))
        _add_to_order(G, order, v, increment)
      end
    end
  
  end

function _create_order(G::AbstractGKM_graph, starting_elem::Int64, rev::Bool)::Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}}
  
  _order::Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}} = Dict{Int64, Tuple{Int64, Int64}}()

  _add_to_order(G, _order, starting_elem, (rev ? -1 : 1))

  return _order
end

function _get_bruhat_order_of_generalized_flag(G::AbstractGKM_graph, rev::Bool = false)::BruhatOrder

  starting_elem::Int64 = rev ? reduce((x,y) -> count('s', G.labels[x]) >= count('s', G.labels[y]) ? x : y, 1:n_vertices(G.g)) : 1

  _order::Dict{Tuple{Int64, Int64}, Vector{Tuple{Int64, Int64}}} = _create_order(G, starting_elem, rev)

  return BruhatOrder(_order, rev, G.labels)

end


@doc raw"""
    get_bruhat_order_of_generalized_flag(R::RootSystem, S::Vector{RootSpaceElem}; descending::Bool=true) -> BruhatOrder

It returns the Bruhat order of the generalized flag variety given by the root system ``R`` with the subset of simple roots given by ``S``. See [`generalized_gkm_flag`](@ref GKMtest.generalized_gkm_flag).
If `descending` is `true`, the Bruhat order is given from the elements of maximal length to the smallest.

# Examples
```jldoctest
julia> R = root_system(:A, 3);

julia> S = simple_roots(R);

julia> get_bruhat_order_of_generalized_flag(R, S[1:2])
Bruhat order in descending order:
Length: 3
  s1*s2*s3 => ["s2*s3"]
Length: 2
  s2*s3 => ["s3"]
Length: 1
  s3 => ["id"]
Length: 0
  id

```
"""
function get_bruhat_order_of_generalized_flag(R::RootSystem, S::Vector{RootSpaceElem}; descending::Bool=true)::BruhatOrder

  base = generalized_gkm_flag(R, S)
  
  return _get_bruhat_order_of_generalized_flag(base, descending)
end


@doc raw"""
    get_bruhat_order_of_generalized_flag(R::RootSystem, indices_of_S; descending::Bool=true) -> BruhatOrder

Same as before, but indicating the indices of the roots in ``S`` instead of the roots itself.

# Examples
```jldoctest
julia> R = root_system(:A, 3);

julia> get_bruhat_order_of_generalized_flag(R, [1]; descending = false)
Bruhat order in ascending order:
Length: 0
  id => ["s2", "s3"]
Length: 1
  s2 => ["s1*s2", "s3*s2", "s2*s3"]
  s3 => ["s3*s2", "s2*s3"]
Length: 2
  s2*s3 => ["s2*s3*s2", "s1*s2*s3"]
  s1*s2 => ["s1*s3*s2", "s1*s2*s3"]
  s3*s2 => ["s1*s3*s2", "s2*s3*s2"]
Length: 3
  s1*s2*s3 => ["s1*s2*s3*s2"]
  s2*s3*s2 => ["s2*s1*s3*s2", "s1*s2*s3*s2"]
  s1*s3*s2 => ["s2*s1*s3*s2", "s1*s2*s3*s2"]
Length: 4
  s2*s1*s3*s2 => ["s1*s2*s1*s3*s2"]
  s1*s2*s3*s2 => ["s1*s2*s1*s3*s2"]
Length: 5
  s1*s2*s1*s3*s2

```
"""
function get_bruhat_order_of_generalized_flag(R::RootSystem, indices_of_S::Union{UnitRange{Int64}, Vector{Int64}}=Int64[]; descending::Bool=true)::BruhatOrder

  base = generalized_gkm_flag(R, indices_of_S)
    
  return _get_bruhat_order_of_generalized_flag(base, descending)
end


function Base.show(io::IO, BO::BruhatOrder)

  if Oscar.is_terse(io)
      # no nested printing
    print(io, "Bruhat order")
  else
      # nested printing allowed, preferably terse
      print(io, "Bruhat order in "*(BO.reversed ? "descending" : "ascending" )*" order")
  end
end
  
  # detailed show
function Base.show(io::IO, ::MIME"text/plain", BO::BruhatOrder)
  
  print(io, "Bruhat order in "*(BO.reversed ? "descending" : "ascending" )*" order:")

  max_length = maximum(lab -> lab[2], keys(BO.order))

  seq = BO.reversed ? (max_length:-1:0) : (0:max_length)
  
  for ord in seq
    print(io, "\nLength: $ord")
    for lab in keys(BO.order)
      if lab[2] == ord
        if ord == last(seq)
          print(io, "\n  $(BO.labels[lab[1]])")
        else
          print(io, "\n  $(BO.labels[lab[1]]) => $([BO.labels[a[1]] for a in BO.order[lab]])")
        end
      end
    end
  end
  
end

@doc raw"""
    generalized_gkm_schubert(R::RootSystem, indices_of_S::Vector{RootSpaceElem}, pt::String) -> AbstractGKM_subgraph

Let ``G`` be the generalized flag variety given by the root system ``R`` with the subset of simple roots given by ``S``. See [`generalized_gkm_flag`](@ref GKMtest.generalized_gkm_flag).
This functions returns the subgraph of the variety ``G`` given by all Schubert cells corresponding to the points less or equal to `pt` in the Bruhat order.

# Examples
```jldoctest generalized_gkm_schubert
julia> R = root_system(:B, 2);

julia> generalized_gkm_schubert(R, "s1*s2")
GKM subgraph of:
GKM graph with 8 nodes, valency 4 and axial function:
s1 -> id => (-1, 1)
s2*s1 -> s1 => (-1, 0)
s1*s2*s1 -> id => (-1, 0)
s1*s2*s1 -> s2*s1 => (-1, -1)
s2 -> id => (0, -1)
s2 -> s2*s1 => (1, -1)
s1*s2 -> s1 => (0, -1)
s1*s2 -> s1*s2*s1 => (1, -1)
s1*s2 -> s2 => (-1, -1)
s2*s1*s2 -> id => (-1, -1)
s2*s1*s2 -> s2*s1 => (0, -1)
s2*s1*s2 -> s1*s2 => (-1, 0)
s1*s2*s1*s2 -> s1 => (-1, -1)
s1*s2*s1*s2 -> s1*s2*s1 => (0, -1)
s1*s2*s1*s2 -> s2 => (-1, 0)
s1*s2*s1*s2 -> s2*s1*s2 => (-1, 1)
Subgraph:
GKM graph with 4 nodes, valency 2 and axial function:
s1 -> id => (-1, 1)
s2 -> id => (0, -1)
s1*s2 -> s1 => (0, -1)
s1*s2 -> s2 => (-1, -1)
```

As before, the subset S can be an subset of simple roots or a subset of indices.
```jldoctest generalized_gkm_schubert
julia> generalized_gkm_schubert(R, [1], "s2")
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
s2 -> id => (0, -1)
s1*s2 -> id => (-1, 0)
s1*s2 -> s2 => (-1, -1)
s2*s1*s2 -> id => (-1, -1)
s2*s1*s2 -> s2 => (-1, 0)
s2*s1*s2 -> s1*s2 => (-1, 0)
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
s2 -> id => (0, -1)

julia> S = simple_roots(R);

julia> generalized_gkm_schubert(R, [S[2]], "s1")
GKM subgraph of:
GKM graph with 4 nodes, valency 3 and axial function:
s1 -> id => (-1, 1)
s2*s1 -> id => (-1, -1)
s2*s1 -> s1 => (-1, 0)
s1*s2*s1 -> id => (-1, 0)
s1*s2*s1 -> s1 => (-1, -1)
s1*s2*s1 -> s2*s1 => (-1, -1)
Subgraph:
GKM graph with 2 nodes, valency 1 and axial function:
s1 -> id => (-1, 1)

```
"""
function generalized_gkm_schubert(R::RootSystem, S::Vector{RootSpaceElem}, pt::String)::AbstractGKM_subgraph

  base = generalized_gkm_flag(R, S)
  @req pt in base.labels "$pt is not a valid label"

  return _generalized_gkm_schubert(base, pt, true)
end

@doc raw"""
    generalized_gkm_schubert(R::RootSystem, indices_of_S::Union{UnitRange{Int64}, Vector{Int64}}, pt::String)  -> AbstractGKM_subgraph
"""
function generalized_gkm_schubert(R::RootSystem, indices_of_S::Union{UnitRange{Int64}, Vector{Int64}}, pt::String)::AbstractGKM_subgraph

  base = generalized_gkm_flag(R, indices_of_S)
  @req pt in base.labels "$pt is not a valid label"
  
  return _generalized_gkm_schubert(base, pt, true)
end

@doc raw"""
    generalized_gkm_schubert(R::RootSystem, pt::String) -> AbstractGKM_subgraph
"""
function generalized_gkm_schubert(R::RootSystem, pt::String)::AbstractGKM_subgraph

    base = generalized_gkm_flag(R, Int64[])
    @req pt in base.labels "$pt is not a valid label"
    
    return _generalized_gkm_schubert(base, pt, true)
  end

function _generalized_gkm_schubert(base::AbstractGKM_graph, pt::String, rev::Bool=true)::AbstractGKM_subgraph

  starting_elem = findfirst(v -> base.labels[v] == pt, 1:n_vertices(base.g))
  suborder = _create_order(base, starting_elem, rev)
  

  vertices::Vector{Int64} = Int64[starting_elem]
  for lab in keys(suborder)
    vertices = vcat(vertices, [a[1] for a in suborder[lab]])
  end

  return gkm_subgraph_from_vertices(base, vertices)
end
