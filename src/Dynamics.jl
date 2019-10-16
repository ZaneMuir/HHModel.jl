# module Dynamics
abstract type AbstractIonChannel end #TODO: generalize simpel and complex channel into generic channels.
abstract type AbstractKinetics end #TODO: extend to more kinetic parameters

@doc raw"""
## Kinetics of an ion channel variable

$m^n$

activation and inactivation function would be encoded
as Boltzmann function:

$$m_\infty(V) = \frac{1 - \zeta}{1 + \exp(\frac{V_{\text{half}} - V}{k})} + \zeta$$

the time constant would be encoded as guassian function:

$$\tau(V) = C_{\text{base}} + C_{\text{amp}} \exp(\frac{- (V_{\text{max}} - V) ^ 2}{\sigma ^ 2})$$

---

### Kinetics(n, Vhalf, k, Cbase, Camp, Vmax, sigma; zeta, state)
- regular variable as conductance model describes.

### Kinetics(n, Vhalf, k; zeta, state, _tau, _type)
- for variable with constant tau
- or instantaneous variable, which $m = m_\infty$ (with `_type = :instantaneous`)

### Kinetics(n, infty, tau)
- custom infty and tau functions

### Kinetics()
- returns an empty kinetic object
"""
mutable struct Kinetics <: AbstractKinetics
    _type::Symbol # :evolving, :intantaneous, :empty
    _state::Symbol # :activation, :inactivation, :custom
    n::Integer

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
    Kinetics(n::Integer, Vhalf::Real, k::Real, Cbase::Real, Camp::Real, Vmax::Real, sigma::Real; zeta::Real=0, state = :activation) = begin
        _sign = state == :activation ? 1 : -1
        _infty = (V) -> (1 - zeta) / (1 + exp(_sign * (Vhalf - V) / k)) + zeta
        _tau = (V) -> Cbase + Camp * exp(- ((Vmax - V) / sigma)^2)
        new(:evolving, state, n, Vhalf, k, zeta, Cbase, Camp, Vmax, sigma, _infty, _tau)
    end

    "kinetics with constant tau or instantaneous kinetic"
    Kinetics(n::Integer, Vhalf::Real, k::Real; zeta::Real=0, _tau = (x)->1, _type=:evolving, state = :activation) = begin
        _sign = state == :activation ? 1 : -1
        _infty = (V) -> (1 - zeta) / (1 + exp(_sign * (Vhalf - V) / k)) + zeta
        new(_type, state, n, Vhalf, k, zeta, nothing, nothing, nothing, nothing, _infty, _tau)
    end

    "custom infty and tau functions"
    Kinetics(n::Integer, infty::Function, tau::Function) = begin
        new(:evolving, :custom,  n, nothing, nothing, nothing, nothing, nothing, nothing, nothing, infty, tau)
    end

    "empty kinetics"
    Kinetics() = begin
        new(:empty, :custom, 0, nothing, nothing, nothing, nothing, nothing, nothing, nothing, (x)->1, (x)->1)
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

include("./Channel/SimpleChannel.jl")
include("./Channel/ComplexChannel.jl")
include("./Channel/GenericChannel.jl")

@doc raw"""
dof(channels::Vector{AbstractIonChannel})

get degree of freedom from model.
only variables of ion channels are counted.
`V` and `I` not included.
"""
function dof(channels::Vector{T}) where {T <: AbstractIonChannel}
    _ans = 0
    for item in channels
        _ans += dof(item) |> sum
    end
    _ans
end

@doc raw"""
simple conductance model, current clamp:

$$C\frac{dV}{dt} = I_{\text{stimulus}} - (\sum_i I_i) $$

using SimpleIonChannels to model a cell;
return an anonymous function that can be used by DifferentialEquations.jl.
"""
function simpleConductanceModel(channels::Vector{T}, stim::Function; C::Real=1) where {T <: AbstractIonChannel}
    nchannel = length(channels)
    nvar = dof(channels)

    return (du, u, p, t) -> begin
        not = (x) -> !x
        v = u[1]
        param = u[2:end-1]
        var_idx = 1

        _current = zeros(nchannel)
        for (idx, item) in enumerate(channels)
            _var_step = sum(dof(item))
            (_icurrent, _iderivitate) = step(item, V=v, var=param[var_idx:var_idx-1+_var_step], E=p.E[item.ion])
            _current[idx] = _icurrent
            du[1+var_idx:var_idx+_var_step] = _iderivitate[not.(isnothing.(_iderivitate))]
            var_idx += _var_step
        end

        I = stim(t, p.stim)
        du[1] = (I - sum(_current)) / C
        u[end] = I

        du, u, p, t
    end
end

function setup_init(channels::Vector{T}, v0::Real) where {T<:AbstractIonChannel}
    _init_val = zeros(dof(channels)+2)
    _init_val[1] = v0
    var_idx = 2

    for (idx, item) in enumerate(channels)
        for each in itr_kinetics(item)
            if each._type == :evolving
                _init_val[var_idx] = each.infty(v0)
                var_idx += 1
            end
        end
    end
    _init_val
end

function current_decompose(solution::ODESolution, model::Vector{T}, tspan::StepRangeLen, param::NamedTuple) where {T<:AbstractIonChannel}
    var = hcat(solution(tspan).u...)
    var_idx = 2
    result = Dict{String, Vector{Float64}}()
    for ch in model
        _var_step = sum(dof(ch))
        result[ch.name] = current(ch, V=var[1, :], var=var[var_idx:var_idx-1+_var_step, :], E=param.E[ch.ion])
        var_idx += _var_step
    end
    result["voltage"] = var[1, :]
    result
end

# voltage clamp
@doc raw"""
simple conductance model, simple voltage clamp:

$$I_{\text{stimulus}} = (\sum_i I_i) $$

return an anonymous function that can be used by DifferentialEquations.jl.
"""
function simpleVoltageClamp(channels::Vector{T}, stim::Function; C::Real=1) where {T <: AbstractIonChannel}
    nchannel = length(channels)
    nvar = dof(channels)

    return (du, u, p, t) -> begin
        not = (x) -> !x
        v = stim(t, p.stim)
        param = u[2:end-1]
        var_idx = 1

        _current = zeros(nchannel)
        for (idx, item) in enumerate(channels)
            _var_step = sum(dof(item))
            (_icurrent, _iderivitate) = step(item, V=v, var=param[var_idx:var_idx-1+_var_step], E=p.E[item.ion])
            _current[idx] = _icurrent
            du[1+var_idx:var_idx+_var_step] = _iderivitate[not.(isnothing.(_iderivitate))]
            var_idx += _var_step
        end

        u[end] = sum(_current)
        u[1] = v

        du, u, p, t
    end
end

#function run()

# end #module
