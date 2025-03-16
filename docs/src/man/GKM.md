# Generalities on GKM graphs

Here we collect the functions that allows the creation of GKM graphs.

# Constructors
These functions allow the construction of GKM varieties.

```@docs
gkm_graph
empty_gkm_graph
GKMadd_edge!
convert_weights
GKM_isValid
```

These are the main properties of GKM graphs.
```@docs
bettiNumbers
valency
rank_torus
is2_indep
is3_indep
```

# Operators
```@docs
GKMsubgraph_from_vertices
blowupGKM
*
```
# Famous examples
These functions allow the construction of GKM varieties.
```@docs
generalized_gkm_flag
flag_variety
grassmannian
gkm_graph_of_toric
```

# Connections
Functions regarding connections.
```@docs
get_GKM_connection
set_GKM_connection!
GKM_isValidConnection
```

# Index

```@index
```