function [w_isentropic, n_stage, cooling_duty, comp_ratio_adj] = compressor_station(p_in, p_out, comp_ratio, gas, T)
% [w_isentropic, n_stage, cooling_duty, comp_ratio_adj] = compressor_station(p_in, p_out, comp_ratio, gas, T)
% returns values in J/mol for isentropic work and cooling duty
% All input values are in SI
% the rest is self explanatory
n_stage=round(log(p_out/p_in)/log(comp_ratio)); % number of compressors
comp_ratio_adj = exp(log(p_out/p_in)/n_stage);  % adjusted compression ratio accourding to the rounded number of compressors
w_isentropic = 0;
cooling_duty = 0;
for i = 1:n_stage
    p1 = p_in*comp_ratio_adj^(i-1);
    p2 = p_in*comp_ratio_adj^i;
    S_in=py.CoolProp.CoolProp.PropsSI("SMOLAR", "T", T, "P", p1, gas); % molar entropy J/(mol.K)
    H_in=py.CoolProp.CoolProp.PropsSI("HMOLAR", "T", T, "P", p1, gas); % molar enthalpy J/mol
    try
        T_out=py.CoolProp.CoolProp.PropsSI("T", "SMOLAR", S_in, "P", p2, gas);
    catch
        f = @(T)(py.CoolProp.CoolProp.PropsSI("SMOLAR", "T", T, "P", p2, gas)-S_in);
        T_out = fzero(f, T);
    end
    H_out=py.CoolProp.CoolProp.PropsSI("HMOLAR", "T", T_out, "P", p2, gas); % molar enthalpy J/mol
    w_isentropic = w_isentropic + H_out-H_in; % J/mol
    
    cooling_duty = cooling_duty+H_out-py.CoolProp.CoolProp.PropsSI("HMOLAR", "T", T, "P", p2, gas); % J/mol
end

end