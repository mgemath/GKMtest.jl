export isrationally_smooth, isregular_word, issmooth_schubert_at_v, issmooth_schubert
function isrationally_smooth(base::AbstractGKM_graph, pt::String, rev::Bool=true)
    
  starting_elem = findfirst(v -> base.labels[v] == pt, 1:n_vertices(base.g))
  suborder = _create_order(base, starting_elem, rev, true)

  expected_valency = length(suborder[starting_elem, count('s', base.labels[starting_elem])])

  (expected_valency == count('s', pt)) || return false
#   println(suborder)

  for lab in keys(suborder)

    # number inbound
    inbound = count(t -> Base.in(lab, suborder[t]), keys(suborder))

    #number outbound
    outbound = length(suborder[lab])

    (inbound + outbound) != expected_valency && return false
    # println("$lab, $inbound $outbound")
  end

  return true
end

# Lemma 5.7
function isregular_word(w::WeylGroupElem)
    pos_roots = positive_roots(root_system(parent(w)))

    l_w = length(w)
    for x in parent(w)
        if x < w

            card_set = count(alpha -> (x < reflection(alpha)*x <= w), pos_roots)

            (card_set > l_w - length(x)) && return false
        end
    end

    return true
end

function issmooth_schubert(w::WeylGroupElem)

  return issmooth_schubert_at_v(w, one(parent(w)))
end
# 7.2.1 Theorem in Singular Loci of Schubert Varieties, Sara Billey and V. Lakshmibai
function issmooth_schubert_at_v(w::WeylGroupElem, v::WeylGroupElem)
  
  @req parent(w) == parent(v) "Weyl groups mismatch"

  Phi = root_system(parent(w))
  s = reflection.(simple_roots(Phi))
  p = length(w)
  word_w = word(w)

  R, a = polynomial_ring(ZZ, :a => (1:rank(Phi)))
  S = fraction_field(R)

  cwv = zero(S) # 7.1.5 Theorem.

  for expon in Iterators.product([0:1 for _ in 1:p]...)
    
    # check if expon is acceptable
    prod(i -> s[word_w[i]]^expon[i], 1:p; init = one(parent(w))) != v && continue
    
    # start the sum
    summand = one(R)
    
    for i in 1:p
      to_add = simple_root(Phi, Int64(word_w[i])) * prod(k -> s[word_w[k]]^expon[k], i:-1:1)
      coef = matrix(ZZ, coefficients(to_add))
      summand *= sum(l -> coef[l]*a[l[2]], eachindex(coef))
    end

    cwv += inv(S(summand))

  end

  cwv *= (-1)^p


  right_pol = ((-1)^(p - length(v))) * one(S)

  for beta in positive_roots(Phi)
    !(reflection(beta) * v <= w) && continue
    
    coef = matrix(ZZ, coefficients(beta))

    right_pol *= inv(S(sum(l -> coef[l]*a[l[2]], eachindex(coef))))

  end

  return cwv == right_pol
end

# function kostant_kumar_polynomial(w::WeylGroupElem, v::WeylGroupElem)

#   @req parent(w) == parent(v) "Weyl groups mismatch"

#   R, a = polynomial_ring(ZZ, :a => (1:rank(root_system(parent(w)))))

#   return _kostant_kumar_polynomial(w, v, R, a, one(R))

# end

# function _kostant_kumar_polynomial(_w::WeylGroupElem, _v::WeylGroupElem, _R::ZZMPolyRing, _a::Vector{ZZMPolyRingElem}, prev_ans::ZZMPolyRingElem)

#   if !(_w <= _v)
#     return zero(_R)
#   end

#   Weyl = parent(_w)
#   Phi = root_system(Weyl)

#   ans::ZZMPolyRingElem = prev_ans
  
#   if _w == _v
#     inversion_set = intersect(positive_roots(Phi), negative_roots(Phi) .* _w)
    
#     for i in 1:rank(Phi)
#       if positive_root(Phi, i) in inversion_set
#         ans *= _a[i]
#       end
#     end
#     return ans
#   end


#   # look for any simple reflection s such that s*w < w
#   s = one(Weyl)
#   alpha = positive_root(Phi, 1)
  
#   for _alpha in positive_roots(Phi)
#     if reflection(_alpha) * _w < _w
#       alpha = _alpha
#       s = reflection(alpha)
#       break
#     end
#   end

#   KK = _kostant_kumar_polynomial(s*_w, _v, _R, _a, ans); println(KK)

#   if _v < s*_v
#     return KK
#   end

#   # here we have v > s*v
#   coef = Int64.(coefficients(alpha*(s*_w)))

#   first_p = sum(i -> _a[coef[i]], eachindex(coef))
  
#   return KK + first_p*_kostant_kumar_polynomial(s*_w, s*_v, _R, _a, ans)

# end

# function kostant_kumar_polynomial(W::WeylGroup, v::WeylGroupElem, w::WeylGroupElem)
#   # Step 1: Check if v ≤ w in Bruhat order
#   if !leq(v, w)  # v not less than or equal to w
#       return 0
#   end

#   # Step 2: Get the list of simple reflections
#   # S = simple_reflections(W)
#   S = gens(W)
#   l_w = length(w)
#   l_v = length(v)

#   # Step 3: Find a reduced expression for w
#   # reduced_w = reduced_word(w)
#   reduced_w = word(w)

#   # Step 4: Compute the polynomial using divided difference operators
#   # Initialize to 1
#   p = 1

#   # Apply divided difference operators according to the steps from w to v
#   for i in 1:l_w
#       s = S[reduced_w[i]]
#       ws = apply(s, w)
#       if leq(v, ws)
#           # divided difference operator: δ_s(p) = (p - s(p)) / α_s
#           alpha_s = simple_root(W, s)
#           p = (p - s(p)) / alpha_s
#           w = ws
#       end
#   end

#   return p
# end
