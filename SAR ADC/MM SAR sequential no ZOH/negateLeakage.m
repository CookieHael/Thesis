function coeff = negateLeakage(sensing_matrix, sensingOrderMatrix, gain, leakage_current, hold_cap, clk_period, coeff)

    nb_MAC = size(sensing_matrix,2);
    leak_cycles_array = zeros(1,length(coeff));
    for i=1:length(coeff)
        idx = find(sensing_matrix(i, :));
        idx2 = find(sensingOrderMatrix(i,:));
        if(isempty(idx))
            idx=384;
        end
        if (isempty(idx2))
            idx2=384;
        end
        leak_cycles = idx2(1)-idx(1);
        leak_cycles_array(i) = leak_cycles;
        coeff(i) = coeff(i) + 0.85*gain*leakage_current*clk_period*leak_cycles/hold_cap;
    end
% coeff = coeff + gain*leakage_current*clk_period*nb_MAC/hold_cap;
fprintf("Average leakage cycles " + num2str(round(mean(leak_cycles_array))) + " cycles .\n");
end
