#-----------------------------------------------------------------------------------------------------
# Spin Operators
#-----------------------------------------------------------------------------------------------------
# Dictionary
function spin_dict(D::Integer)
    J = (D-1)/2
    coeff = [sqrt(i*(D-i)) for i = 1:D-1]
    sp = sparse(1:D-1, 2:D, coeff, D, D)
    sm = sp'
    sx = (sp+sm)/2
    sy = (sp-sm)/2
    sz = sparse(1:D, 1:D, J:-1:-J)
    s0 = sparse(1.0I, D, D)
    Dict([('+', sp), ('-', sm), ('x', sx), ('y', sy), ('z', sz), ('1', s0)])
end
# Atomic spin matrix
function spin_atom(s::String, dic::Dict)
    n = length(s)
    ny = sum(si == 'y' for si in s)
    temp = n == 1 ? dic[s[1]] : kron([dic[si] for si in s]...)
    P = mod(ny, 2) == 0 ? (-1)^(ny÷2) : (1im)^ny
    P * temp
end
# General spin matrix
function spin(spins; D::Integer=2)
    dic = spin_dict(D)
    res = sum(ci * spin_atom(si, dic) for (ci, si) in spins)
    Array(res)
end

#---------------------------------------------------------------------------------------------------
# Other properties
#---------------------------------------------------------------------------------------------------
export entropy
entropy(S::AbstractVector; cutoff::Real=1e-14) = -sum(i*log(i) for i in S if i > cutoff)

export inner_product
function inner_product(T1::Array, T2::Array)
    α, β = size(T1, 3), size(T2, 3)
    K = Kraus(T1, conj(T2), :r)
    ρ = rand(ComplexF64, α, β) |> vec
    e, v = eigsolve(K, ρ)
    abs(e[1])
end

inner_product(T::Array) = inner_product(T, T)
