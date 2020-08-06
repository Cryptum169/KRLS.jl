eye(n) = 1.0 * Matrix(I, n ,n);

export FilterParam, KRLSFilter

struct FilterParam
    order::Int
    fALD::Float64
    iDtMax::Int
end

mutable struct KRLSFilter
    mKinv::Array
    vAlpha::Array
    mP::Array
    mDict::Array
    m::Int
    mode::Char
    params::FilterParam
    function KRLSFilter(vx::Array, vy::Array, mode::Char, params::FilterParam)
        # Frst Check vx and vy shape, should be a vertical vector
        ktt = vx * vx'
        mKinv = 1 ./ ktt
        vAlpha = vx * vy'
        # mP = eye(length(vx))
        mP = I
        mDict = vx
        m = length(vx)
        new(mKinv, vAlpha, mP, mDict, m, mode, params)
    end
end

@doc raw"""
	update(vx::Array, vy::Array, filter::KRLSFilter)

Update Filter intrinsic parameters for new observations.

# Examples
```julia-repl
julia> update(vx, vy, filter)
```
"""
function update(vx::Array, vy::Array, filter::KRLSFilter)
    at = filter.mKinv * kernel(vx, filter.m, filter.mode, filter.mDict);
    ktt = [filter.mDict vx]' * vx;
    deltat = abs(ktt[end] - dot(kernel(vx, filter.m, filter.mode, filter.mDict), at))
    # Dot is required here to ensure it produces a scalar

    rowC = size(filter.mDict, 2);

    if deltat > filter.params.fALD && rowC < filter.params.iDtMax
        filter.mDict = [filter.mDict vx];
        filter.mKinv = 1 / deltat * [deltat * filter.mKinv + at * at' -at; -at' 1];
        (r_mp, c_mp) = size(filter.mP);
        filter.mP = [filter.mP zeros(r_mp, 1); zeros(1, c_mp) 1];

        filter.vAlpha = [filter.vAlpha + at / deltat * (vy' - kernel(vx, filter.m, filter.mode,
            filter.mDict)' * filter.vAlpha); 1 / deltat * (vy' - kernel(vx, filter.m, filter.mode,
            filter.mDict)' * filter.vAlpha)];
        filter.m = filter.m + 1;
    else
        qt = filter.mP * at / (1 + (at' * filter.mP * at)[1]);
        filter.mP = filter.mP - qt * at' * filter.mP;
        filter.vAlpha = filter.vAlpha + filter.mKinv * qt * (vy' - kernel(vx, filter.m, filter.mode,
            filter.mDict)' * filter.vAlpha);
    end

    return filter;
end

function predict(vx::Array, filter::KRLSFilter)
    (kernel(vx, filter.m, filter.mode, filter.mDict)' * filter.vAlpha)';
end
