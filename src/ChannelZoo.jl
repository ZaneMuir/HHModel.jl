@doc raw"""
## low voltage gated potassiumd channel

### simple version

$$i_{\text{ltk}} = g_{\text{ltk}} w^4 z (V - E_{\text{K}})$$
- subtype:
    - `kv1`: $\zeta=0.5$
    - `ikcnq`: $\zeta=1.0$
    - `ia`: $\zeta=0.0$

<!-- TODO: ### complex version -->

- modified from `inf_tau_w_ltk_rm.m`, `inf_tau_z_ltk_rm.m` and `I_ltk_rm.m` 

"""
function low_voltage_gated_potassium(g::Real; iklcomplex=false, subtype=:kv1)
    _param = (
        kv1 =   (0.5, 1000), 
        ikcnq = (1.0, 10000),
        ia =    (0.0, 10)
    )
    _zeta, _camp = _param[subtype]
    
    if iklcomplex
        #TODO
        return SimpleIonChannel("low-voltage-gated potassium, complex", :potassium, 
            g, 0, 0,
            Kinetics(), Kinetics())
    else
        # Change back to 100 for kv1 like
        _w_tau = (V) -> 100*(6*exp((V+60)/6)+16*exp(-(V+60)/45))^(-1)+1.5 # IKL 
        _w = Kinetics(-44.5, 8.4, 0, _tau = _w_tau)
        
        _z_tau = (V) -> _camp*(exp((V+60)/20)+exp(-(V+60)/8))^(-1)+50;
        _z = Kinetics(-71.0, 10.0, _zeta, _tau = _z_tau, state=:inactivation)

        return SimpleIonChannel("low-voltage-gated potassium, simple", :potassium, 
            g, 4, 1,
            _w, _z)
    end
end

@doc raw"""
## high voltage gated potassium channel

"""
function high_voltage_gated_potassium(g::Real)
    # TODO
end

@doc raw"""
## hudgkin huxley sodium current

$$i_{\text{Na}} = g_{\text{Na}} m^3 h (V - E_{\text{Na}})$$

- modified from `inf_tau_m_rm.m`, `inf_tau_h_rm.m` and `I_na_rm.m`
"""
function hh_sodium(g::Real)
    _m_tau = (V) -> 10 / (5*exp((V+60)/18)+36*exp(-(V+60)/25))+0.04
    _m = Kinetics(-38.0, 7.0, 0, _tau = _m_tau)
    
    _h_tau = (V) -> 100 / (7*exp((V+60)/11)+10*exp(-(V+60)/25))+0.6
    _h = Kinetics(-65.0, 6.0, 0, _tau=_h_tau, state=:inactivation)
    
    SimpleIonChannel("hh sodium", :sodium, 
        g, 3, 1,
        _m, _h)
end

@doc raw"""
## hh potassium channel

$$i_{\text{K}} = g_{\text{K}} a ^ 4 b c (V - E_{\text{K}})$$
"""
function hh_potassium(g::Real)
    #TODO
end

@doc raw"""
## Ih current
$$i_{\text{h}} = g_{\text{h}} r (V - E_{\text{h}})$$

- modified from `inf_tau_r_rm.m` and `I_h_rm.m`
"""
function ihcurrent(g::Real)
    a=83;
    b = 7.6;
    c = 0.0002;
    d = 5.3;
    e=313;
    _r_tau = (V) -> (1/c) * 1/(exp(-(V+a)/b) + exp((V+a)/d))+e;
    _r = Kinetics(-100.0, 7.0, 0, _tau=_r_tau, state=:inactivation)
    
    SimpleIonChannel("ih current", :ih, 
        g, 0, 1,
        Kinetics(), _r)
end

@doc raw"""
## leakage current
$$i_{\text{leak}} = g_{\text{leak}} (V - E_{\text{leak}})$$
"""
function leakage(g::Real)
    SimpleIonChannel("leakage", :leak, 
        g, 0, 0,
        Kinetics(), Kinetics())
end