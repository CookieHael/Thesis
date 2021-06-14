function [SRBM] = generateSRBM(M, N, S)
%GENERATESRBM Generate MxN (1,s)-sparse random binary matrix

SRBM = zeros(M,N);
for i=1:N
    ones = randperm(M, S);
    SRBM(ones,i) = 1;
end
r = rank(SRBM);
    
while(r ~= M && M>=N && S~=M)
    SRBM = zeros(M,N);
    for i=1:N
        ones = randperm(M, S);
        SRBM(ones,i) = 1;
    end
    r = rank(SRBM);
end

end

