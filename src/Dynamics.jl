# module Dynamics
abstract type AbstractIonChannel end

@doc raw"""
## Kinetics of an ion channel variable

activation and inactivation function would be encoded
as Boltzmann function:

$$m_\infty(V) = \frac{1 - \zeta}{1 + \exp(\frac{V_{\text{half}} - V}{k})} + \zeta$$

the time constant would be encoded as guassian function:

$$\tau(V) = C_{\text{base}} + C_{\text{amp}} \exp(\frac{- (V_{\text{max}} - V) ^ 2}{\sigma ^ 2})$$

---

### Kinetics(Vhalf, k, Cbase, Camp, Vmax, sigma)
- regular variable as conductance model describes.

### Kinetics(Vhalf::Real, k::Real; _tau = (x)->1, _type=:evolving)
- for variable with constant tau
- or instantaneous variable, which $m = m_\infty$ (with `_type = :instantaneous`)

### Kinetics()
- returns an empty kinetic object
"""
mutable struct Kinetics
    _type::Symbol #:evolving, :intantaneous, :empty
    
    Vhalf::Union{Nothing, Real}
    k::Union{Nothing, Real}
    zeta::Union{Nothing, Real}
    
    Cbase::Union{Nothing, Real}
    Camp::Union{Nothing, Real}
    Vmax::Union{Nothing, Real}
    sigma::Union{Nothing, Real}
    
    infty::Function
    tau::Function
    
    "set kinetics"
    Kinetics(Vhalf::Real, k::Real, zeta::Real, Cbase::Real, Camp::Real, Vmax::Real, sigma::Real; state = :activation) = begin
        _sign = state == :activation ? 1 : -1
        _infty = (V) -> (1 - zeta) / (1 + exp(_sign * (Vhalf - V) / k)) + zeta
        _tau = (V) -> Cbase + Camp * exp(- ((Vmax - V) / sigma)^2)
        new(:evolving, Vhalf, k, zeta, Cbase, Camp, Vmax, sigma, _infty, _tau)
    end
    
    "kinetics with constant tau or instantaneous kinetic"
    Kinetics(Vhalf::Real, k::Real, zeta::Real; _tau = (x)->1, _type=:evolving, state = :activation) = begin
        _sign = state == :activation ? 1 : -1
        _infty = (V) -> (1 - zeta) / (1 + exp(_sign * (Vhalf - V) / k)) + zeta
        new(_type, Vhalf, k, zeta, nothing, nothing, nothing, nothing, _infty, _tau)
    end
    
    "custom infty and tau functions"
    Kinetics(infty::Function, tau::Function) = begin
        new(:evolving, nothing, nothing, nothing, nothing, nothing, nothing, nothing, infty, tau)
    end
    
    "empty kinetics"
    Kinetics() = begin
        new(:empty, nothing, nothing, nothing, nothing, nothing, nothing, nothing, (x)->1, (x)->1)
    end
end

function activation(knt::Kinetics, V::Real)
    if isnothing(knt.zeta)
        return knt.infty(V)
    else
        return (1 - knt.zeta) / (1 + exp((knt.Vhalf - V) / knt.k)) + knt.zeta
    end
end

function inactivation(knt::Kinetics, V::Real)
    if isnothing(knt.zeta)
        return knt.infty(V)
    else
        return (1 - knt.zeta) / (1 + exp(- (knt.Vhalf - V) / knt.k)) + knt.zeta
    end
end

function time_constant(knt::Kinetics, V::Real)
    if isnothing(knt.Cbase)
        return knt.tau(V)
    else
        return knt.Cbase + knt.Camp * exp(- ((knt.Vmax - V) / knt.sigma)^2)
    end
end

@doc raw"""
conductance based channel

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
    a::Integer
    b::Integer
    m::Kinetics
    h::Kinetics
end

@doc raw"""
"""
mutable struct ComplexIonChannel <: AbstractIonChannel
    #TODO
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

@doc raw"""
conductance based channel

$$i = g m^a h^b (V - E)$$
$$\dot{m} = \frac{m_\infty(V) - m}{\tau_m(V)}$$
$$\dot{h} = \frac{h_\infty(V) - h}{\tau_h(V)}$$
"""
function step(ch::SimpleIonChannel; V::Real, m::Real, h::Real, E::Real)
    i = ch.g * (V - E) 
    
    if ch.m._type == :evolving
        i = i * m ^ ch.a
        dm = (ch.m.infty(V) - m) / ch.m.tau(V)
    elseif ch.m._type == :instantaneous
        i = i * ch.m.infty(V) ^ ch.a
        dm = nothing
    else
        dm = nothing
    end
    
    if ch.h._type == :evolving
        i = i * h ^ ch.b
        dm = (ch.h.infty(V) - h) / ch.h.tau(V)
    elseif ch.h._type == :instantaneous
        i = i * ch.h.infty(V) ^ ch.b
        dh = nothing
    else
        dh = nothing
    end
    
    return (i, dm, dh)
end

@doc raw"""
dof(channels::Vector{SimpleIonChannel})

get degree of freedom from model.
only variables of ion channels are counted.
`V` and `I` not included.
"""
function dof(channels::Vector{SimpleIonChannel})
    counter = 0
    for item in channels
        if item.m._type == :evolving
            counter = counter + 1
        end
        if item.h._type == :evolving
            counter = counter + 1
        end
    end
    counter
end

function dof(channel::SimpleIonChannel)
    _m = channel.m._type == :evolving ? 1 : 0
    _h = channel.h._type == :evolving ? 1 : 0
    (_m, _h)
end

@doc raw"""
simple conductance model:

using SimpleIonChannels to model a cell;
return an anonymous function that can be used by DifferentialEquations.jl.
"""
function simpleConductanceModel(channels::Vector{SimpleIonChannel}, stim::Function)
    nchannel = length(channels)
    nvar = dof(channels)
    
    return (du, u, p, t) -> begin
        v = u[1]
        param = u[2:end-1]
        var_idx = 1
        dvar_idx = 1
        
        _current = zeros(nchannel)
        for (idx, item) in enumerate(channels)
            (_m, _h) = dof(item)
            
            # retrive activation or inactivation variable
            if _m == 1
                _m_val = param[var_idx]
                var_idx = var_idx + 1
            else
                _m_val = 1
            end
            
            if _h == 1
                _h_val = param[var_idx]
                var_idx = var_idx + 1
            else
                _h_val = 1
            end
            
            # step update
            (_item_i, _item_dm, _item_dh) = step(item, V=v, m=_m_val, h=_h_val, E=p.E[item.ion])
            _current[idx] = _item_i
            
            # udpate du value for activation or inactivation variable
            if _m == 1
                du[1 + dvar_idx] = _item_dm
                dvar_idx = dvar_idx + 1
            end
            
            if _h == 1
                du[1 + dvar_idx] = _item_dh
                dvar_idx = dvar_idx + 1
            end
        end
        
        # udpate dV and current input
        du[1] = stim(t) - sum(_current)
        u[end] = stim(t)
        du
    end
end

# end #module