module HHModelCompoenents

function inf_tau_h_rm(V)
    i = (1+exp((V+65)/6))^(-1);
    t = 100*(7*exp((V+60)/11)+10*exp(-(V+60)/25))^(-1)+0.6;
    i, t
end

function inf_tau_m_rm(V)
    i = (1+exp(-(V+38)/7))^(-1);
    t = 10*(5*exp((V+60)/18)+36*exp(-(V+60)/25))^(-1)+0.04;
    i, t
end

function inf_tau_n_rm(V)
    a, b = alpha_beta_n(V);
    i = a / (a+b);
    t = 1 / (a+b);
    i, t
end

function alpha_beta_n(V)
    if V != 60
        a = -0.01*(V+60)/(exp(-(V+60)/10)-1);
    else
        a = 0.1;
    end
    b = 0.125*exp(-(V+70)/80);
    a, b
end

end