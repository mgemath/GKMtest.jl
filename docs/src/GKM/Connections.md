# Connections

Following [GZ98](@cite), a *connection* on a GKM graph $G$ is a function $\Delta_e:E_p \rightarrow E_q$ for each edge $e: p\rightarrow q$ where $E_p$ is the set of oriented edges in $G$ starting at $p$.
It is required to satisfy the following conditions.
 * For each edge $e$, $\Delta_e(e) = \overline{e}$, i.e. the same edge with reverse orientation.
 * If $w:E\rightarrow M$ is the axial function, then for every edge $e:p\rightarrow q$ and $e'\in E_p$, there exists $a\in\mathbb{Z}$ such that $w(e') - w(\Delta_e(e')) = a\cdot w(e)$.

If $G$ is the GKM graph of a GKM variety $X$, then these integers $a$ are the degrees of the equivariant line bundles into which $TX$ splits when restricted to the $\mathbb{P}^1$ represented by $e$.

## Existence and uniqueness of connections

Given a GKM graph $G$ that comes from a GKM variety $X$, it always has a connection for the geometric reason sketched above.
However, it is often convenient if one does not have to specify $\Delta$ manually, so the package does it automatically whenever $\Delta$ is unique by one of the following two reasons.
Sufficient conditions for uniqueness of the connection (if it exists) are:
 * The valency of $G$ is at least 3 and $G$ is $3$-independent, i.e. the weights of every three edges starting at the same vertex are linearly independent.
 * The valency of $G$ is at most 2.

If neither of these two conditions hold and $G$ is not the output of a standard construction, a choice of connection can be specified manually using `set_connection!`.

```@docs
get_GKM_connection
set_GKM_connection!
GKM_isValidConnection
```