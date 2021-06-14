function coeff = negateLeakage(sensing_matrix, sensingOrderMatrix, gain, leakage_current, hold_cap, clk_period, coeff)

     nb_MAC = size(sensing_matrix,2);

    for i=1:length(coeff)
        idx = find(sensing_matrix(i, :));
        idx2 = find(sensingOrderMatrix(i,:));
        if(~isempty(idx))
            if (isempty(idx2))
                idx2=384;
            end
            leak_cycles = idx2(1)-idx(1)-1;
            coeff(i) = coeff(i) + gain*leakage_current*clk_period*leak_cycles/hold_cap;
        end
    end

% coeff = coeff + gain*leakage_current*clk_period*nb_MAC/hold_cap;

end
