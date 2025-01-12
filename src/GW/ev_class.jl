function ev(j::Int64, cc)::EquivariantClass

  rule = :(_ev(dt, $j, $cc))
  return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function _ev(dt::GW_decorated_tree, j::Int64, cc)

  v = imageOf(dt.marks[j], dt)
  
  return cc[v]
end