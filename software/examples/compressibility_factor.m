function Z = compressibility_factor(p, T, comp)
    % function Z = compressibility_factor(p, T, comp)
    % p: pressure [Pa]
    % T: temperature [K]
    % comp: componet name in CoolProp; string type
    % Z: compressibility factor [-]
    % example:
    %  Z = compressibility_factor(1e5, 300, 'Air'); % returns Z for air
    % example components: Air, CH4, water, NH3, ...
    
    Z = py.CoolProp.CoolProp.PropsSI('Z','P',p,'T',T,comp); % [-]
end