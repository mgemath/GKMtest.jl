export planarcurves

function planarcurves(n::Int64, d::Int64, P::EquivariantClass)
  return planarcurves(n, d, [P])
end

function planarcurves(n::Int64, d::Int64, P_input::Array{EquivariantClass})

  G::AbstractGKM_graph = projective_space(AbstractGKM_graph, n)
  beta::CurveClass_type = d*curve_class(G, Edge(1, 2));
  show_bar::Bool = true
  check_degrees::Bool = false

  inputLength = length(P_input)
  inputSize = size(P_input)
  inputKeys = keys(P_input)
  @req inputLength > 0 "gromov_witten needs at least one input for P_input..."

  @req beta != zero(parent(beta)) "Beta must be non-zero"

  H2 = GKM_second_homology(G)
  R = G.equivariantCohomology
  res = zeros(R.coeffRingLocalized, inputSize...)
  if !is_effective(H2, beta)
    return res
  end

  P = [P_input[k].func for k in inputKeys]
  con = get_any_connection(G)
  @req !isnothing(con) "GKM graph needs a connection!"

  ########
  # this part is needed for the generation of colorings
  nc::Dict{Int64,Vector{Int64}} = Dict{Int64,Vector{Int64}}()
  for v in 1:n_vertices(G.g)
    nc[v] = sort(all_neighbors(G.g, v))
  end
  #########

  #########
  # Dict in order to store H
  h_dict::Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}} = Dict{Tuple{Int64, Int64, Int64}, AbstractAlgebra.Generic.FracFieldElem{QQMPolyRingElem}}() # Lambda_gamma_e_dict
  ########

  max_n_vert::Int64 = _max_n_edges(H2, beta) + 1

  if show_bar #set up progress data
    number_trees = A000055(max_n_vert)
    # Count the number of trees with at most max_n_vert vertices and a graph homomorphism to
    # G.g, also counting ways to distribute the marked points among the vertices.
    # This does not cound the edge multiplicities.
    threshold = sum(vert -> number_trees[vert] * n_vertices(G.g) * ((length(nc[1]))^(vert - 1)) * binomial(vert+n_marks-1, n_marks), 2:max_n_vert)
    progress_bar::Progress = Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green)
    current_graph = 0
  end

  # n_marks = length(classes)
  # iterate undecorated trees:
  for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert]) # generation of level sequences
    tree = LStoGraph(ls) # from level sequence to graph
    tree_aut = count_iso(ls)

    CI, parents, subgraph_ends = col_it_init(ls, nc) # generation of colorings
    # iterate maps from tree to G.g:
    for col in CI   # colorings Iterator
      top_aut::Int64 = count_iso(ls, col)

      Multi = _multiplicities(H2, [Edge(col[src(e)], col[dst(e)]) for e in edges(tree)], beta)

      # iterate location of marks on the tree
      for m_inv in Combinatorics.with_replacement_combinations(1:nv(tree), n_marks)
        
        aut = count_iso(ls, col, m_inv)

        # iterate edge multiplicities
        for edgeMult_array in Multi

          PROD = prod(edgeMult_array)
          euler = zero(R.coeffRing)

          edgeMult = Dict{Edge, Int}(edges(tree) .=> edgeMult_array)
          
          # Iterate numbering of the marks on the tree, picking only one per isomorphism class
          # Details here have to do with the colors iterator from Colors.jl.
          for m in Base.Iterators.filter(mul_per -> top_aut == 1 || isempty(mul_per) || maximum(mul_per) < 3 || ismin(ls, col, mul_per, parents, subgraph_ends), multiset_permutations(m_inv, n_marks))

            
            dt = decoratedTree(G, tree, col, edgeMult, m)
            
            Class = [Base.invokelatest(P[k], dt) for k in inputKeys]
            #println("Tree($(dt.vDict)), Marks($(dt.marks)), Mult($(dt.edgeMult)):")
            #println(Class)
            all(c -> c==0, Class) && continue

            if euler == zero(R.coeffRing)
              euler = Euler_inv(dt; check_degree=check_degrees)//(PROD * aut)
              for e in edges(tree)
                triple = (edgeMult[e], min(col[src(e)], col[dst(e)]), max(col[src(e)], col[dst(e)]))
                if !haskey(h_dict, triple)
                    h_dict[triple] = _h(Edge(col[src(e)], col[dst(e)]), triple[1], con, R; check=false, check_degrees=check_degrees)
                end
                euler *= h_dict[triple]
              end
            end
            
            ## Planar Curves part
            
            ##
            res += Class.*euler
          end
        end

        if show_bar #update the progress bar
          current_graph += tree_aut ÷ top_aut
          update!(progress_bar, current_graph,
              showvalues=[(:"Total number of graphs", threshold), (:"Current graph", current_graph)])
        end
      end
    end
  end
  return res
end

# function 