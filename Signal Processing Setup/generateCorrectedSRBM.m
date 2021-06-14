function [SRBM, CorrectedSRBM] = generateCorrectedSRBM(M, N, S, K)
%GENERATESRBM Generate MxN (1,s)-sparse random binary matrix

if (M==S)
    error("Number of ones per column equal to number of columns. Can't generate full-rank matrix.");
elseif (M>N)
    warning("Matrix dimensions incorrect for compression.");
end
rng shuffle
SRBM = zeros(M,N);
for i = 1 : N
    o = randperm(M);
    indices = o(1:S);
    col = zeros(M,1);
    col(indices) = ones(S,1);
    SRBM(:,i) = col;
end  

nb_ones_column = sum(SRBM,2);
CorrectedSRBM = SRBM;
for j=1:M
    nb_ones = nb_ones_column(j);
    for i=1:N
        if (CorrectedSRBM(j,i) == 1)
            CorrectedSRBM(j,i) = (K/(K+1))^(nb_ones-1);
            nb_ones = nb_ones-1;
        end
    end
end

