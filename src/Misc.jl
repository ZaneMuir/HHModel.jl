function detect_cross_pnt(arr, thr, way=:up, gap=1)
    _idx_repo = []
    try
        idx = findall(arr .> thr)
        idx_diff = findall(diff(idx) .> 1)
        _idx_repo = [idx[1]; idx[idx_diff]; idx[idx_diff .+ 1]; idx[end]]
    catch BoundError
        return nothing
    end
    
    if way == :up
        _check = (x, i) -> x[i-1] < thr < x[i] #< x[i+1]
    elseif way == :down
        _check = (x, i) -> x[i-1] > x[i] > x[i+1]
    else
        return nothing
    end
    
    sort!(_idx_repo)
    _result = Vector{Int}()
    _previous = -9999
    
    for idx in _idx_repo
            if _check(arr, idx) && (idx - _previous > gap)
                _result = [_result; idx]
                _previous = idx
            end
    end
    
    _result
end