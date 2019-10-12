@doc raw"""
conductance based channel, generic format

$$i = g\cdot \text{krule}(\text{var}) \cdot (V - E)$$

"""
mutable struct GenericIonChannel <: AbstractIonChannel #TODO: optimize the structure for better and easier use.
    name::String
    ion::Symbol
    g::Real
    
    krule::Function
    kvars::Vector{Kinetics}
end

function dof(ch::GenericIonChannel)
    _ans = zeros(Int, length(ch.kvars))
    for (idx, item) in enumerate(ch.kvars)
        if item._type == :evolving
            _ans[idx] = 1
        end
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

#NOTE: not tested
function current(ch::GenericIonChannel; V::Vector{T}, var::Array{T, 2}, E::T) where {T<:Real}
    _var_idx = 1
    _krule = (var_item) -> ch.krule(ch.kvars, var_item)
    _var_reform = [Tuple(var[:, idx]) for idx =1:size(var, 2)]
    current = ch.g .* _krule.(_var_reform) .* (V .- E)
    return current
end