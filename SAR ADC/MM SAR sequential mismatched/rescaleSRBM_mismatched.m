function CorrectedSRBM = rescaleSRBM_mismatched(SRBM, C1_array, C2_array)

SRBM = SRBM~=0;

CorrectedSRBM = double(SRBM);
for j=1:size(SRBM,1)
    nb_ones = sum(SRBM(j,:));
    K = C2_array(j)/(C2_array(j)+C1_array);
    for i=1:size(SRBM,2)
        if (CorrectedSRBM(j,i) == 1)
            CorrectedSRBM(j,i) = K^(nb_ones-1);
            nb_ones = nb_ones-1;
        end
    end
end



end

