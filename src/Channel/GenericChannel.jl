@doc raw"""
conductance based channel, generic format

$$i = g\cdot \text{krule}(\text{var}) \cdot (V - E)$$

"""
mutable struct GenericIonChannel <: AbstractIonChannel
    name::String
    ion::Symbol
    g::Real
    
    krule::Function
    kvars::Vector{Kinetics}
end

function dof(ch::GenericIonChannel)
    _ans = 0
    for item in ch.kvars
        _ans += item._type == :evolving ? 1 : 0
    end
    _ans
end

function update!(ch::GenericIonChannel)
    for item in ch.kvars
        if ~isnothing(item.Vhalf)
            if item._state == :activation
                item.infty = (V) -> (1 - item.zeta) / (1 + exp((item.Vhalf - V) / item.k)) + item.zeta
            elseif item._state == :inactivation
                item.infty = (V) -> (1 - item.zeta) / (1 + exp(-(item.Vhalf - V) / item.k)) + item.zeta
            else
                nothing
            end
        end

        if ~isnothing(item.Cbase)
            item.tau = (V) -> item.Cbase + item.Camp * exp(- ((item.Vmax - V) / item.sigma)^2)
        end
    end
    ch
end

function step(ch::GenericIonChannel; V::Real, var::Vector{T}, E::Real) where {T<:Real}
    _var_idx = 1
    derivative = zeros(size(var))
    
    current = ch.g * ch.krule(ch.kvars, var) * (V - E)
    for item in ch.kvars
        if item._type == :evolving
            derivative[_var_idx] = (item.infty(V) - var[_var_idx]) / item.tau(V)
            _var_idx += 1
        else
            nothing
        end
    end
    
    return (current, derivative)
end

function itr_kinetics(ch::GenericIonChannel)
    ch.kvars
end