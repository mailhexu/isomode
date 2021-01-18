using SparseArrays
struct ts
 i :: Number
 A:: SparseMatrixCSC
 B:: Array{<:Any, 2}
end

ts(M::Array{<:Any, 2})=ts(1, M, M)

M=[1.0 0;0 1]
println(M)
t=ts(M)
println(typeof(t.A))
println(t.A)
println(typeof(t.B))
println(t.B)
println(t.i)
println(typeof(t.i))
