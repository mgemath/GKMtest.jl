import Base: *, //, +, -, one, zero, ==

_CohomType = Union{AbstractAlgebra.Generic.FreeModuleElem, FreeModElem, Array}

function QH_class(G::AbstractGKM_graph, class)
  H2 = GKM_second_homology(G)
  beta0 = zero(H2.H2)
  nv = n_vertices(G.g)
  coeffs = Dict{CurveClass_type, QH_coeff_type}()
  c0 = zero(G.equivariantCohomology.cohomRingLocalized)
  Rloc = G.equivariantCohomology.coeffRingLocalized
  g = gens(G.equivariantCohomology.cohomRingLocalized)
  for i in 1:nv
    c0 += Rloc(class[i]) * g[i]
  end
  coeffs[beta0] = c0
  return QHRingElem(G, coeffs)
end

function _QH_remove_zero_coeffs!(c::QHRingElem)::QHRingElem
  for b in keys(c.coeffs)
    if is_zero(c.coeffs[b])
      delete!(c.coeffs, b)
    end
  end
  return c
end

function +(c1::QHRingElem, c2::QHRingElem)::QHRingElem
  @req c1.gkm == c2.gkm "QH ring elements don't belong to the same GKM graph"
  res = QHRingElem(c1.gkm, Dict{CurveClass_type, QH_coeff_type}())
  for b in keys(c1.coeffs)
    res.coeffs[b] = c1.coeffs[b]
  end
  for b in keys(c2.coeffs)
    if haskey(res.coeffs, b)
      res.coeffs[b] += c2.coeffs[b]
    else
      res.coeffs[b] = c2.coeffs[b]
    end
  end
  return _QH_remove_zero_coeffs!(res)
end

function +(c1::_CohomType, c2::QHRingElem)::QHRingElem
  c1QH = QH_class(c2.gkm, c1)
  return c1QH + c2
end

function +(c1::QHRingElem, c2::_CohomType)::QHRingElem
  c2QH = QH_class(c1.gkm, c2)
  return c1 + c2QH
end

function +(c1::QHRingElem, c2)::QHRingElem
  return c1 + one(c2) * c2
end

function +(c1, c2::QHRingElem)::QHRingElem
  return one(c1)*c1 + c2
end

# Here, s should be a scalar, i.e. anything that embeds into the localized coefficient ring H_T(*)_loc.
function *(s, c::QHRingElem)::QHRingElem
  res = QHRingElem(c.gkm, Dict{CurveClass_type, QH_coeff_type}())
  for b in keys(c.coeffs)
    res.coeffs[b] = (s//1) * c.coeffs[b]
  end
  return _QH_remove_zero_coeffs!(res)
end

function *(c::QHRingElem, s)::QHRingElem
  return *(s, c)
end

function *(c1::QHRingElem, c2::QHRingElem)::QHRingElem
  @req c1.gkm == c2.gkm "QH ring elements don't belong to the same GKM graph"
  #calculate QH structure constants in all beta. This throws an error if G is not strictly NEF.
  SC = QH_structure_constants(c1.gkm)
  res = QHRingElem(c1.gkm, Dict{CurveClass_type, QH_coeff_type}())
  for b1 in keys(c1.coeffs)
    for b2 in keys(c2.coeffs)
      for b in keys(SC)
        tmp = quantumProduct(c1.gkm, b, c1.coeffs[b1], c2.coeffs[b2])
        k = b1+b2+b
        if haskey(res.coeffs, k)
          res.coeffs[k] += tmp
        else
          res.coeffs[k] = tmp
        end
      end
    end
  end
  return _QH_remove_zero_coeffs!(res)
end

function *(c1::_CohomType, c2::QHRingElem)::QHRingElem
  c1QH = QH_class(c2.gkm, c1)
  return c1QH * c2
end

function *(c1::QHRingElem, c2::_CohomType)::QHRingElem
  c2QH = QH_class(c1.gkm, c2)
  return c1 * c2QH
end

function -(c1::QHRingElem, c2::QHRingElem)::QHRingElem
  @req c1.gkm == c2.gkm "QH ring elements don't belong to the same GKM graph"
  res = QHRingElem(c1.gkm, Dict{CurveClass_type, QH_coeff_type}())
  for b in keys(c1.coeffs)
    res.coeffs[b] = c1.coeffs[b]
  end
  for b in keys(c2.coeffs)
    if haskey(res.coeffs, b)
      res.coeffs[b] -= c2.coeffs[b]
    else
      res.coeffs[b] = -c2.coeffs[b]
    end
  end
  return _QH_remove_zero_coeffs!(res)
end

function -(c1::_CohomType, c2::QHRingElem)::QHRingElem
  c1QH = QH_class(c2.gkm, c1)
  return c1QH - c2
end

function -(c1::QHRingElem, c2::_CohomType)::QHRingElem
  c2QH = QH_class(c1.gkm, c2)
  return c1 - c2QH
end

function -(c1::QHRingElem, c2)::QHRingElem
  return c1 - one(c1) * c2
end

function -(c1, c2::QHRingElem)::QHRingElem
  return one(c2)*c1 - c2
end

function //(c::QHRingElem, s)::QHRingElem
  @req s != 0 "Cannot divide QHRingElem by zero!"
  res = QHRingElem(c.gkm, Dict{CurveClass_type, QH_coeff_type}())
  for b in keys(c.coeffs)
    res.coeffs[b] = (1//s) * c.coeffs[b]
  end
  return res
end

function one(c::QHRingElem)::QHRingElem
  res = zero(c)
  for g in gens(c.gkm.equivariantCohomology.cohomRingLocalized)
    res += g
  end
  return res
end

function zero(c::QHRingElem)::QHRingElem
  res = QHRingElem(c.gkm, Dict{CurveClass_type, QH_coeff_type}())
  res.coeffs[zero(GKM_second_homology(c.gkm).H2)] = zero(c.gkm.equivariantCohomology.cohomRingLocalized)
  return res
end

function ==(c1::QHRingElem, c2::QHRingElem)::Bool
  c1.gkm != c2.gkm && return false
  for b in keys(c1.coeffs)
    if haskey(c2.coeffs, b)
      c1.coeffs[b] != c2.coeffs[b] && return false
    else
      c1.coeffs[b] != 0 && return false
    end
  end
  for b in keys(c2.coeffs)
    if haskey(c1.coeffs, b)
      c2.coeffs[b] != c1.coeffs[b] && return false
    else
      c2.coeffs[b] != 0 && return false
    end
  end
  return true
end

function ==(c1::QHRingElem, c2::_CohomType)::Bool
  c2QH = QH_class(c1.gkm, c2)
  return c1 == c2QH
end

function ==(c1::_CohomType, c2::QHRingElem)::Bool
  return c2 == c1
end

function ==(c1::QHRingElem, c2)::Bool
  return c1 == one(c1)*c2
end

function ==(c1, c2::QHRingElem)::Bool
  return c2 == c1
end

function Base.show(io::IO, c::QHRingElem)

  if isempty(c.coeffs)
    print(io, "0")
    return
  end

  l = length(keys(c.coeffs))
  counter = 0
  # no nested printing
  for b in keys(c.coeffs)
    counter += 1
    print(io, c.coeffs[b])
    print(io, " q^")
    print(io, b)
    if counter != l
      if Oscar.is_terse(io)
        print(io, " + ")
      else
        print(io, "\n + ")
      end
    end
  end
end

# detailed show
function Base.show(io::IO, ::MIME"text/plain", c::QHRingElem)
  show(io, c)
end