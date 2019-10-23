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

export SimpleIonChannel, GenericIonChannel, ComplexIonChannel, update!, dof
export VoltageClampSimulation, CurrentClampSimulation
export low_voltage_gated_potassium, high_voltage_gated_potassium, hh_sodium, hh_potassium, ihcurrent, leakage
export export_as_atf, export_as_hdf5, export_as_mat


end
