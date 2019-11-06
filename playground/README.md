# Demos

## Get Started
if you are using the `Julia` for the first time.
Run the `setup.jl` script to install dependent pacakges.

```bash
> cd playground
> julia ./setup.jl
```

## Contents
- SimplifiedModel.ipynb
    - modeling simplified "Persistent Sodium Plus Potassium Model"
        - parameters from Izhikevich book
    - setting up biophysical parameters
    - setting up stimulus function
    - running a current clamp simulation
    - visualizing the voltage and current trace
- HudgkinHuxleyModel.ipynb
    - modeling standard Hudgkin Huxley Model
        - parameters from Izhikevich book
    - setting up biophysical parameters and stimulus function
    - running a voltage clamp simulation
    - visualizing the membrane voltage and current trace
- SingleCompartment_RM.ipynb
    - modeling full model described in _Hight and Kalluri, 2016_.
        - channels are built in `ChannelZoo.jl`
    - setting up biophysical parameters and stimulus function
    - running a current clamp protocol with steps
    - visualizing the voltage and current trace