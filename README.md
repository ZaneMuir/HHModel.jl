# HHModel.jl
![](http://www.wtfpl.net/wp-content/uploads/2012/12/wtfpl-badge-4.png)

conductance based neuronal modeling.


## demos
![SingleCompartmentDemo](./playground/SingleCompartmentDemo.ipynb)

## references
- [Hight, A. E. & Kalluri, R. A biophysical model examining the role of low-voltage-activated potassium currents in shaping the responses of vestibular ganglion neurons. J. Neurophysiol. 116, 503–521 (2016).](https://doi.org/10.1152/jn.00107.2016)

## TODO list
- documentation and workflow
- complex ltk channel
- post-synaptic potential/current simulation
- model struct and wrap simulation into one method.

## change logs

#### 10/15/2019
- make playground files more readable.
- verification with Maltab codes
- bug fix:
    - currents of channels
    - initial values

#### 10/11/2019
- voltage clamp mode
- calculate current of each channel during current clamp

#### 10/07/2019
- basic framework
    - Kinetics
    - SimpleIonChannel
    - ComplexIonChannel
    - GenericIonChannel
    - modelgeneration
- channel zoo
    - high voltage gated potassium channel
    - low voltage gated potassium channel
    - hh sodium channel
    - hh potassium channel
    - ih current
    - leakage current
