function mu = viscosity(p, T, comp)
    % function mu = viscosity(p, T, comp)
    % return the viscosity of a fluid at a constant temperature and  pressure
    % p: pressure [Pa]
    % T: temperature [K]
    % comp: componet name in CoolProp; string type
    % mu: viscosity [Pa.s]
    % comp is a string that contains the name of a CoolProp fluid
    mu = py.CoolProp.CoolProp.PropsSI('VISCOSITY','P',p,'T',T,comp); % Pa.s
end
