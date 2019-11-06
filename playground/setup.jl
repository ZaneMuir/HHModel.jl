using Pkg

Pkg.add.([
        "DifferentialEquations",
        "Plots", # using GR. Or install other backends by yourself.
        "MAT",
        "HDF5",
        "IJulia",
        "MultivariateStats", # for PCA analysis
        "ProgressMeter"
        ])

using IJulia # install IJulia kernel
