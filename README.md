# HHModel.jl
![](http://www.wtfpl.net/wp-content/uploads/2012/12/wtfpl-badge-4.png)

conductance based neuronal modeling, using [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl)

## demos
check out jupyter notebooks in `playground` folder.

## references

- ChannelZool.jl: [Hight, A. E. & Kalluri, R. A biophysical model examining the role of low-voltage-activated potassium currents in shaping the responses of vestibular ganglion neurons. J. Neurophysiol. 116, 503–521 (2016).](https://doi.org/10.1152/jn.00107.2016)
- basic framework: [Izhikevich, E. Dynamical Systems in Neuoscience: The Geometry of Excitability and Bursting. The MIT Press, 2007](http://www.izhikevich.org/publications/dsn/index.htm)

## TODO list
- documentation and workflow
- complex ltk channel
- model struct and wrap simulation into one method, or through Ensemble Simulations
- intrinsically bursting neuronal model
- chattering neuronal model
- neural network

