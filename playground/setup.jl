using Pkg

Pkg.add.(["DifferentialEquations",
        "Plots",
        "MAT",
        "HDF5",
        "IJulia",
        "MultivariateStats",
        "ProgressMeter",
        "FFTW",
        "DSP"])

using IJulia # install IJulia kernel
