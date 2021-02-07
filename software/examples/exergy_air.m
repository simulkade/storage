function ex = exergy(p, T, comp)
    % ex = exergy_air(p, T)
    % p: pressure [Pa]
    % T: temperature [K]
    % returns exergy of comp in J/kg
    T0 = 15+273.15; K
    p0 = 1e5; % Pa

    H = py.CoolProp.CoolProp.PropsSI('H','P',p,'T',T,comp); % J/mol
    S = py.CoolProp.CoolProp.PropsSI('S','P',p,'T',T, comp); % J/mol.K
    H0 = py.CoolProp.CoolProp.PropsSI('H','P',p0,'T',T0,comp); % J/mol
    S0 = py.CoolProp.CoolProp.PropsSI('S','P',p0,'T',T0,comp); % J/mol.K
    ex = (H-H0)-T0*(S-S0);
end