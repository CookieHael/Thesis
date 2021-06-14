function SRBM = correctMatrix(SRBM, K)
%CONVERTMATRIX Summary of this function goes here
%   Detailed explanation goes here
nb_ones_column = sum(SRBM,2);
for j=1:size(SRBM,1)
    nb_ones = nb_ones_column(j);
    for i=1:size(SRBM,2)
        if (SRBM(j,i) == 1)
            SRBM(j,i) = (K/(K+1))^(nb_ones-1);
            nb_ones = nb_ones-1;
        end
    end
end

