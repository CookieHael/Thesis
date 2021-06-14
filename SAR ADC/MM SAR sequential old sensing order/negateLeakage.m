function coeff = negateLeakage(sensing_matrix, gain, leakage_current, hold_cap, clk_period, coeff)

     nb_MAC = size(sensing_matrix,2);

    for i=1:length(coeff)
        idx = find(sensing_matrix(i, :));
        if(~isempty(idx))
            leak_cycles = nb_MAC-idx(1)+1;
            coeff(i) = coeff(i) + gain/2*leakage_current*clk_period*leak_cycles/hold_cap;
        end
    end

% coeff = coeff + gain*leakage_current*clk_period*nb_MAC/hold_cap;

end
