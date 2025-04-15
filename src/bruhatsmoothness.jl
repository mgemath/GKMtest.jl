export isrationally_smooth, isregular_word
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
