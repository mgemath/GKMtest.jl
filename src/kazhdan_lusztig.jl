export R_polynomial, kazhdan_lusztig

function R_polynomial(v::WeylGroupElem, w::WeylGroupElem)
  @req parent(w) == parent(v) "Weyl groups mismatch"

  A, q = polynomial_ring(ZZ, "q")

  return _R_polynomial(v, w, A, q)
end

function _R_polynomial(v::WeylGroupElem, w::WeylGroupElem, A::ZZPolyRing, q::ZZPolyRingElem)

  w == v && return one(A)
  
  v < w || return zero(A)

  for s in gens(parent(w))
    if w*s < w
      if v*s < v
        return _R_polynomial(v*s, w*s, A, q)
      end
      if v < v*s
        return (q-one(A))*_R_polynomial(v, w*s, A, q) + q*_R_polynomial(v*s, w*s, A, q)
      end
    end
  end
    
end

function kazhdan_lusztig(v::WeylGroupElem, w::WeylGroupElem)
  @req parent(w) == parent(v) "Weyl groups mismatch"

  A, q = polynomial_ring(ZZ, "q")

  return _kazhdan_lusztig(v, w, A, q)
end

function _kazhdan_lusztig(v::WeylGroupElem, w::WeylGroupElem, A::ZZPolyRing, q::ZZPolyRingElem)
  
  v <= w || return zero(A)

  deg_P = div(length(w) - length(v) - 1, 2)

  (deg_P < 1) && return one(A)

  right_hand_side = zero(A)

  for y in parent(w)
    if v < y && y <= w
      right_hand_side += ((_R_polynomial(v, y, A, q) * _kazhdan_lusztig(y, w, A, q)) % q^(deg_P+1))  # cut terms of degree >deg_P
    end
  end

  return -right_hand_side
end
