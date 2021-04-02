function CorrectedSRBM = rescaleSRBM(SRBM, K)

SRBM = SRBM~=0;

nb_ones_column = sum(SRBM,2);
CorrectedSRBM = SRBM;
for j=1:size(SRBM,1)
    nb_ones = nb_ones_column(j);
    for i=1:size(SRBM,2)
        if (CorrectedSRBM(j,i) == 1)
            CorrectedSRBM(j,i) = (K/(K+1))^(nb_ones-1);
            nb_ones = nb_ones-1;
        end
    end
end

end

