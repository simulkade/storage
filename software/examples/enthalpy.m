function rho = exergy_air(p, T, comp)
    % function rho = density(p, T, comp)
    % p: pressure [Pa]
    % T: temperature [K]
    % comp: componet name in CoolProp; string type
    % rho: mass density [kg/m3]
    % example:
    %  rho_water = density(1e5, 300, 'water'); % returns density of water
    % example components: Air, CH4, water, NH3, ...

    rho = py.CoolProp.CoolProp.PropsSI('D','P',p,'T',T,comp); % kg/m3
end