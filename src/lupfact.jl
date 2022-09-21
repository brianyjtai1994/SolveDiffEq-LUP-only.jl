# LUP: Factorization
function lupf!(lu::MatIO, p::VecIO{Int}) # @code_warntype ✓
    nrow, ncol = size(lu)
    @inbounds p[end] = ncol
    for kx in 1:ncol # column-wise
        #### Find pivoting row
        kpiv = kx  # kth column pivot
        amax = 0.0 # max. abs. of kth column
        for ix in kx:nrow
            @inbounds temp = abs(lu[ix,kx])
            if temp > amax
                kpiv = ix
                amax = temp
            end
        end
        @inbounds p[kx] = kpiv
        #### Check singularity
        iszero(amax) && error("lupf!: Singular matrix.")
        #### Row-swap
        if kx ≠ kpiv
            for jx in eachindex(1:ncol)
                swap!(lu, kx, jx, kpiv, jx) # interchange
            end
            @inbounds p[end] += 1
        end
        #### Scale the 1st column
        @inbounds lukkinv = inv(lu[kx,kx])
        iszero(lukkinv) && (lukkinv = eps())
        @simd for ix in kx+1:nrow
            @inbounds lu[ix,kx] *= lukkinv
        end
        #### Update the rest block
        @inbounds for jx in kx+1:ncol, ix in kx+1:nrow
            lu[ix,jx] -= lu[ix,kx] * lu[kx,jx]
        end
    end
    return nothing
end

# `trsv!`: upper-triangular, no-transposed, normal-diagoanl
function trsv!(::Val{:U}, ::Val{:N}, ::Val{:N}, A::MatI, x::VecIO)
    n = length(x)
    n ≠ size(A, 2) && error("trsv!: n ≠ size(A, 2)")
    iszero(n) && return nothing
    for j in n:-1:1
        if !iszero(x[j])
            @inbounds x[j] /= A[j,j]
            @inbounds temp = x[j]
            @simd for i in j-1:-1:1
                @inbounds x[i] -= temp * A[i,j]
            end
        end
    end
    return nothing
end

# `trsv!`: lower-triangular, no-transposed, unit-diagoanl
function trsv!(::Val{:L}, ::Val{:N}, ::Val{:U}, A::MatI, x::VecIO)
    n = length(x)
    n ≠ size(A, 2) && error("trsv!: n ≠ size(A, 2)")
    iszero(n) && return nothing
    for j in 1:n
        if !iszero(x[j])
            @inbounds temp = x[j]
            @simd for i in j+1:n
                @inbounds x[i] -= temp * A[i,j]
            end
        end
    end
    return nothing
end

# LUP: Solve Linear System
function lups!(lu::MatI, x::VecIO, p::VecI{Int})
    @inbounds for i in eachindex(x)
        swap!(x, i, p[i])
    end
    trsv!(VALL, VALN, VALU, lu, x)
    trsv!(VALU, VALN, VALN, lu, x)
    return nothing
end
