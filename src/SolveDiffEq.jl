module SolveDiffEq

const VecI  = AbstractVector
const VecIO = AbstractVector
const MatI  = AbstractMatrix
const MatIO = AbstractMatrix

const VALL = Val(:L)
const VALN = Val(:N)
const VALU = Val(:U)

# @code_warntype ✓
function swap!(v::VecIO, i::Int, j::Int)
    @inbounds temp = v[i]
    @inbounds v[i] = v[j]
    @inbounds v[j] = temp
    return nothing
end

# @code_warntype ✓
function swap!(m::MatIO, i1::Int, j1::Int, i2::Int, j2::Int)
    @inbounds temp     = m[i1,j1]
    @inbounds m[i1,j1] = m[i2,j2]
    @inbounds m[i2,j2] = temp
    return nothing
end

include("./lupfact.jl")

end # module SolveDiffEq
