function class_one()::EquivariantClass

    rule = :(_class_one(dt))

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function _class_one(dt::GW_decorated_tree)

    return 1
end