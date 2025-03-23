# Generalities on GKM graphs

## Definition

GKM varities have been introduced in [GKM98](@cite). They are defined as smooth complex varieties with a torus action such that:
* The action has a finite number of fixed points, and a finite number of 1-dimensional orbits.
* The odd cohomology of the varities is trivial.

We associate to a GKM variety $X$ with torus $T$ a GKM graph. This is the datum of:
* A graph having the fixed points as vertices, such that two vertices are connected by an edge if there is a 1-dimensional orbit passing through the two fixed points.
* An axial function $\mathrm{w}\colon E \rightarrow M$ from the set of oriented edges of the graph to the weight lattice of $T$.

In this package, the axial function can take values in a free $\mathbb{Z}$-module or in a free $\mathbb{Q}$-module.

Famous examples of GKM varities are projective space, see [Standard Constructions](STDconstructions.md) for more examples.

## Index

```@index
```