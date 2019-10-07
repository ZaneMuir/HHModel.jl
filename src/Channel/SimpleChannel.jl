@doc raw"""
conductance based channel, simple format

$$i = \bar{g} m^a h^b (V - E)$$

$$i = g m^a h^b (V - E)$$
$$\dot{m} = \frac{m_\infty(V) - m}{\tau_m(V)}$$
$$\dot{h} = \frac{h_\infty(V) - h}{\tau_h(V)}$$

activation and inactivation function would be encoded
as Boltzmann function:

$$m_\infty(V) = \frac{1 - \zeta}{1 + \exp(\frac{V_{\text{half}} - V}{k})} + \zeta$$

the time constant would be encoded as guassian function:

$$\tau(V) = C_{\text{base}} + C_{\text{amp}} \exp(\frac{- (V_{\text{max}} - V) ^ 2}{\sigma ^ 2})$$
"""
mutable struct SimpleIonChannel <: AbstractIonChannel
    name::String
    ion::Symbol
    g::Real
    
    m::Kinetics
    h::Kinetics
end

function update!(ch::SimpleIonChannel)
    if ~isnothing(ch.m.Vhalf)
        ch.m.infty = (V) -> (1 - ch.m.zeta) / (1 + exp((ch.m.Vhalf - V) / ch.m.k)) + ch.m.zeta
    end

    if ~isnothing(ch.m.Cbase)
        ch.m.tau = (V) -> ch.m.Cbase + ch.m.Camp * exp(- ((ch.m.Vmax - V) / ch.m.sigma)^2)
    end

    if ~isnothing(ch.h.Vhalf)
        ch.h.infty = (V) -> (1 - ch.h.zeta) / (1 + exp(-(ch.h.Vhalf - V) / ch.h.k)) + ch.h.zeta
    end

    if ~isnothing(ch.h.Cbase)
        ch.h.tau = (V) -> ch.h.Cbase + ch.h.Camp * exp(- ((ch.h.Vmax - V) / ch.h.sigma)^2)
    end

    ch
end

function dof(channel::SimpleIonChannel)
    _m = channel.m._type == :evolving ? 1 : 0
    _h = channel.h._type == :evolving ? 1 : 0
    [_m, _h]
end

function step(ch::SimpleIonChannel; V::Real, var::Vector{T}, E::Real) where {T <: Real}
    _var_idx = 1
    i = ch.g * (V - E)

    if ch.m._type == :evolving
        m = var[_var_idx]
        _var_idx += 1
        
        i = i * m ^ ch.m.n
        dm = (ch.m.infty(V) - m) / ch.m.tau(V)
    elseif ch.m._type == :instantaneous
        i = i * ch.m.infty(V) ^ ch.m.n
        dm = nothing
    else
        dm = nothing
    end

    if ch.h._type == :evolving
        h = var[_var_idx]
        
        i = i * h ^ ch.h.n
        dh = (ch.h.infty(V) - h) / ch.h.tau(V)
    elseif ch.h._type == :instantaneous
        i = i * ch.h.infty(V) ^ ch.h.n
        dh = nothing
    else
        dh = nothing
    end

    return (i, [dm, dh])
end

function itr_kinetics(ch::SimpleIonChannel)
    [ch.m, ch.h]
end