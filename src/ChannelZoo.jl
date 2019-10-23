@doc raw"""
## low voltage gated potassiumd channel

### simple version

$$i_{\text{ltk}} = g_{\text{ltk}} w^4 z (V - E_{\text{K}})$$
- subtype:
    - `kv1`: $\zeta=0.5$
    - `ikcnq`: $\zeta=1.0$
    - `ia`: $\zeta=0.0$

<!-- TODO: ### complex version -->

#### ref
- Hight and Kalluri, 2016
"""
function low_voltage_gated_potassium(g::Real; iltkcomplex=false, subtype=:kv1)
    _param = (
        kv1 =   (0.5, 1000),
        ikcnq = (1.0, 10000),
        ia =    (0.0, 10)
    )
    _zeta, _camp = _param[subtype]

    if iltkcomplex
        #TODO
        return SimpleIonChannel("ltk_cmplx", :potassium,
            g, 0, 0,
            Kinetics(), Kinetics())
    else
        # Change back to 100 for kv1 like
        _w_tau = (V) -> 100*(6*exp((V+60)/6)+16*exp(-(V+60)/45))^(-1)+1.5 # IKL
        _w_infty = (V) -> (1+exp(-(V+44.5)/8.4))^(-1/4)
        _w = Kinetics(4, _w_infty, _w_tau)

        _z_tau = (V) -> _camp*(exp((V+60)/20)+exp(-(V+60)/8))^(-1)+50;
        _z = Kinetics(1, -71.0, 10.0, zeta=_zeta, _tau=_z_tau, state=:inactivation)

        return SimpleIonChannel("ltk", :potassium,
            g, _w, _z)
    end
end

@doc raw"""
## high voltage gated potassium channel

#### ref
- Hight and Kalluri, 2016
"""
function high_voltage_gated_potassium(g::Real; phi=0.85)
    _n_tau = (V) -> 100*(11*exp((V+60)/24)+21*exp(-(V+60)/23))^(-1)+0.7;
    _n_infty = (V) -> (1+exp(-(V+15)/5))^(-1/2)
    _n = Kinetics(2, _n_infty, _n_tau)

    _p_tau = (V) -> 100*(4*exp((V+60)/32)+5*exp(-(V+60)/22))^(-1)+5;
    _p = Kinetics(1, -23, 6, _tau=_p_tau)

    ComplexIonChannel("htk", :potassium,
           g, [phi, 1-phi],
           [(_n, Kinetics()),
            (_p, Kinetics())]
    )
end

@doc raw"""
## hudgkin huxley sodium current

$$i_{\text{Na}} = g_{\text{Na}} m^3 h (V - E_{\text{Na}})$$

- modified from `inf_tau_m_rm.m`, `inf_tau_h_rm.m` and `I_na_rm.m`

#### ref
- Hight and Kalluri, 2016
"""
function hh_sodium(g::Real)
    _m_tau = (V) -> 10 / (5*exp((V+60)/18)+36*exp(-(V+60)/25))+0.04
    _m = Kinetics(3, -38.0, 7.0, _tau = _m_tau)

    _h_tau = (V) -> 100 / (7*exp((V+60)/11)+10*exp(-(V+60)/25))+0.6
    _h = Kinetics(1, -65.0, 6.0, _tau=_h_tau, state=:inactivation)

    SimpleIonChannel("ina", :sodium,
        g, _m, _h)
end

@doc raw"""
## hh potassium channel

$$i_{\text{K}} = g_{\text{K}} a ^ 4 b c (V - E_{\text{K}})$$

#### ref
- Hight and Kalluri, 2016
"""
function hh_potassium(g::Real)
    ik_rule = (knt, var) -> begin
        (a,b,c) = var
        a ^ 4 * b * c
    end

    _a_tau = (V) -> 100*(7*exp((V+60)/14)+29*exp(-(V+60)/24))^(-1)+0.1
    _a_infty = (V) -> (1+exp(-(V+31)/6))^(-1/4)
    _a = Kinetics(4, _a_infty, _a_tau)

    _b_tau = (V) -> 1000*(14*exp((V+60)/27)+29*exp(-(V+60)/24))^(-1)+1;
    _b_infty = (V) -> (1+exp((V+66)/7))^(-1/2);
    _b = Kinetics(1, _b_infty, _b_tau)

    _c_tau = (V) -> 90*(1+exp(-(V+66)/17))^(-1)+10;
    _c_infty = (V) -> (1+exp((V+66)/7))^(-1/2);
    _c = Kinetics(1, _c_infty, _c_tau)

    GenericIonChannel("ik", :potassium,
        g, ik_rule,
        [_a, _b, _c]
    );
end

@doc raw"""
## Ih current
$$i_{\text{h}} = g_{\text{h}} r (V - E_{\text{h}})$$

#### ref
- Hight and Kalluri, 2016
- for different activation level: Kalluri's Lab
"""
function ihcurrent(g::Real; level=1)
    _activation_level = Dict(
        1 => (86.86183, 5.76769, 3393.34171, 6.67393, 295.95798, -96.11, 8.1),
        2 => (85.3077,6.34179,3183.003,8.72304,274.3381,-91.67525,8.2),
        3 => (83.75372,6.34179,2972.664705,10.77215,252.71829,-89.2405,8.3),
        4 => (83.23572,6.437473,2902.551813,11.456875,246.67834,-88.428416,8.35),
        5 => (82.71772,6.533156,2832.43892,12.1416,240.63839,-87.6163,8.40),
        6 => (82.19972,6.62884,2762.32603,12.826325,234.598445,-86.80425,8.45),
        7 => (80.64572,6.91589,2551.9877,14.8805,209.4786,-84.368,8.6),
    )
    a,b,c,d,e,_vhalf,_k = _activation_level[level]
    # _r_tau = (V) -> (1/c) * 1/(exp(-(V+a)/b) + exp((V+a)/d))+e;
    _r_tau = (V) -> ((1/c)*1/(exp(-(V+a)/b)) + (exp((V+a)/d)))+e;
    _r = Kinetics(1, _vhalf, _k, _tau=_r_tau, state=:inactivation)

    SimpleIonChannel("ih", :ih,
        g, Kinetics(), _r)
end

@doc raw"""
## leakage current
$$i_{\text{leak}} = g_{\text{leak}} (V - E_{\text{leak}})$$

"""
function leakage(g::Real)
    SimpleIonChannel("leakage", :leak,
        g, Kinetics(), Kinetics())
end