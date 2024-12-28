###############################################################################
#
#   GKMtest
#
###############################################################################

# GKM_graph{C<:RingElem} = Tuple{Graph{Directed}, Dict{Edge, AbstractAlgebra.Generic.FreeModule{C}}}

@attributes mutable struct AbstractGKM_graph{C<:RingElem} #<: GKM_graph{C}
  g::Graph{Directed}
  labels::Vector
  α::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{C}}

  function AbstractGKM_graph{C}(
    g::Graph{Directed},
    labels::Vector,
    α::Dict{Edge, AbstractAlgebra.Generic.FreeModuleElem{C}};
    check::Bool=true,
  ) where {C<:RingElement}

    if check
    # perform all checks
    end
    return new{C}(g, labels, α)
  end
end