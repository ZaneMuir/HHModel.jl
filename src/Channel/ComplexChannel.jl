@doc raw"""
conductance based channel, complex format

$$i = g\cdot ( \sum_j \phi m^a h^b ) \cdot (V - E)$$

"""
mutable struct ComplexIonChannel <: AbstractIonChannel
    name::String
    ion::Symbol
    g::Real
    
    weights::Vector{Real}
    var::Vector{Tuple{Kinetics, Kinetics}}
end

function dof(channel::ComplexIonChannel)
    _ans = zeros(Int, 4)
    _idx = 0
    for item in channel.var
        for each in item
            _idx += 1
            if each._type == :evolving
                _ans[_idx] = 1
            end
        end
    end
    _ans
end

function update!(ch::ComplexIonChannel)
    for item in ch.var
        m = item[1]
        if ~isnothing(m.Vhalf)
            m.infty = (V) -> (1 - m.zeta) / (1 + exp((m.Vhalf - V) / m.k)) + m.zeta
        end

        if ~isnothing(m.Cbase)
            m.tau = (V) -> m.Cbase + m.Camp * exp(- ((m.Vmax - V) / m.sigma)^2)
        end
        
        h = item[1]
        if ~isnothing(h.Vhalf)
            h.infty = (V) -> (1 - h.zeta) / (1 + exp(- (h.Vhalf - V) / h.k)) + h.zeta
        end

        if ~isnothing(h.Cbase)
            h.tau = (V) -> h.Cbase + h.Camp * exp(- ((h.Vmax - V) / hsigma)^2)
        end
    end
    ch
end

function step(ch::ComplexIonChannel; V::Real, var::Vector{T}, E::Real) where {T <: Real}
    _var_idx = 1
    derivative = zeros(size(var))
    _kinetics = 0
    
    for (idx, item) in enumerate(ch.var)
        _tmp = 1
        for each in item
            if each._type == :evolving
                _var_item = var[_var_idx]
                

                _tmp *= _var_item ^ each.n
                derivative[_var_idx] = (each.infty(V) - _var_item) / each.tau(V)
                
                _var_idx += 1
            elseif each._type == :instantaneous
                _tmp *= each.infty(V) ^ each.n
            else
                nothing
            end
        end
        _kinetics += ch.weights[idx] * _tmp
    end
    
    current = ch.g * _kinetics * (V - E)
    return (current, derivative)
end

function itr_kinetics(ch::ComplexIonChannel)
    _ans = Vector{Kinetics}(undef, length(ch.var)*2)
    for (idx, item) in enumerate(ch.var)
        _ans[idx*2-1] = item[1]
        _ans[idx*2] = item[2]
    end
    _ans
end

function current(ch::ComplexIonChannel; V::Vector{T}, var::Array{T, 2}, E::T) where {T <: Real}
    _var_idx = 1
    _kinetics = zeros(size(V))
    
    for (idx, item) in enumerate(ch.var)
        _tmp = ones(size(V))
        for each in item
            if each._type == :evolving
                _var_item = var[_var_idx, :]
                _tmp .*= _var_item .^ each.n
                _var_idx += 1
            elseif each._type == :instantaneous
                _tmp .*= each.infty.(V) .^ each.n
            else
                nothing
            end
        end
        _kinetics .+= ch.weights[idx] .* _tmp
    end
    
    current = ch.g .* _kinetics .* (V .- E)
    current
end