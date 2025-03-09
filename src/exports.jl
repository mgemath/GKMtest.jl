# This file contains all exports statements for the GKMtest module.

export GKM_graph #abstract type
export gkm_graph
export GKM_initialize!
export valency
export rank_torus
export connection
export GKMproj_space
export is3_indep
export is2_indep
export GKMadd_edge!
export empty_gkm_graph
# export flag_gkm_graph
export GKM_isValid
export edgeFromLabels

# GKMconnections.jl
export get_GKM_connection
export set_GKM_connection!
export build_GKm_connection
export connection_a_from_con
export connection_map_from_a
export GKM_isValidConnection

# cohomology.jl
export isGKMclass
export weightClass
export scalar, zero, one, multiply, eulerClass, pointClass, integrateClass, integrateGKMClass, PDClass
export firstChernClass

# betti.jl
export bettiNumbers

# GKMsubgraphs.jl
export GKMsubgraph_from_vertices
export GKMsubgraph_from_edges
export GKM_isValidSubgraph
export isCompatible

# blowup.jl
export blowupGKM

# curveClasses.jl
export GKM_second_homology
export edgeCurveClass
export all_classes
export isEffectiveCurveClass
export chernNumber
export isStrictlyNEF