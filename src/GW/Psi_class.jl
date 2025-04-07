function Psi(a)::EquivariantClass

    rule = :(_Psi(dt, $a))
    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function Psi(a::Int...)::EquivariantClass

    rule = :(_Psi(dt, $a))
    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function _Psi(dt::GW_decorated_tree, a::Int64...)
    return _Psi(dt, [a])
end

function _Psi(dt::GW_decorated_tree, a::Tuple{Vararg{Int64}})
    return _Psi(dt, [a...])
end

function _Psi(dt::GW_decorated_tree, a::Vector{Int64})

    findfirst(x -> x > 0, a) === nothing && return 1 # F(1) #if all of them are zero or a is empty
    g = dt.tree
    marks = dt.marks


    ans = 1 #F(1)

    # local q1::fmpq = fmpq(1)
    # local temp1::fmpq = fmpq(1)
    local Sum_ai::Int64
    local n::Int64
    local M::Int64
    # local d = Dict(edges(g) .=> weights) #assign weights to edges
    local inv_marks::Dict{Int64,Vector{Int64}} = invert_marks(marks, nv(g))

    for v in 1:nv(g)

        a_v = Int64[]
        for i in inv_marks[v]
            (i > length(a) || a[i] == 0) && continue
            push!(a_v, a[i])
        end

        Sum_ai = sum(a_v)
        Sum_ai == 0 && continue #if S contains only zeros, or it is empty, continue

        n = length(all_neighbors(g, v)) + length(inv_marks[v])

        n > 2 && Sum_ai > n - 3 && return 0# F(0)

        #If no previous condition holds, then n>1
        if n == 2 #necessary |S_v| == 1
            M = (-1)^a_v[1]
        else # n>2 and Sum_ai <= n - 3
            M = multinomial((n - 3 - Sum_ai, a_v...,))
        end


        local s1 = 0# F(0)

        for w in all_neighbors(g, v)
            e = Edge(v,w)
            # wev = weight_class(imageOf(e, dt), R) // edgeMult(e, dt)
            s1 += edgeMult(e, dt) // weight_class(imageOf(e, dt), dt.gkm) #  1 // wev 
        end
        ans *= M * (s1^(-Sum_ai))
    end

    return ans
end