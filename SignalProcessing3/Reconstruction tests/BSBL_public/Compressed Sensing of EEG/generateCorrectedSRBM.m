function [SRBM] = generateCorrectedSRBM(M, N, S, K)
%GENERATESRBM Generate MxN (1,s)-sparse random binary matrix

if (M==S)
    error("Number of ones per column equal to number of columns. Can't generate full-rank matrix.");
elseif (M>N)
    warning("Matrix dimensions incorrect for compression.");
end

r = 0;
SRBM = zeros(M,N);
while (~(r == M))
    SRBM = zeros(M,N);
    for i=1:N
        ones_indices = randperm(M, S);
        SRBM(ones_indices,i) = 1;
    end
    r = rank(SRBM);
end
nb_ones_column = sum(SRBM,2);
for j=1:M
    nb_ones = nb_ones_column(j);
    for i=1:N
        if (SRBM(j,i) == 1)
            SRBM(j,i) = (K/(K+1))^(nb_ones-1);
            nb_ones = nb_ones-1;
        end
    end
end

