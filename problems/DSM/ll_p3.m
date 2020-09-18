function p3 = ll_p3(xu, xl)
n_lvar = size(xl, 2);
n_uvar = size(xu, 2);
p3 = 30 * abs(sin(pi/n_lvar.* (xl(:, 2:n_lvar) - xu(:, 2:n_uvar))));
end