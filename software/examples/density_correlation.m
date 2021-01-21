function rho_cor = density_correlation(p1, p2, T, comp)
    % function rho_cor = density_correlation(p1, p2, T, comp)
    % returns a density correlation for component comp at constant temperature 
    % and in pressure range [p1, p2]
    p = linspace(p1, p2, 50);
    rho = arrayfun(@(x)density(x, T, comp), p);
    cf = polyfit(p-p1, log(rho/rho(1)), 1);
    cf = cf(1);
    rho_cor = @(p)(rho(1)*exp(cf*(p-p1)));
end