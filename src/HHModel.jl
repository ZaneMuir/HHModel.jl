module HHModel

using Printf
using DifferentialEquations

include("./Dynamics.jl")
include("./ChannelZoo.jl")
include("./Misc.jl")

using MAT
using HDF5
using Dates
include("./IO.jl")

end
