function [SRBM] = generateCorrectedSRBM(M, N, S)
%GENERATESRBM Generate MxN (1,s)-sparse random binary matrix
r = 0;
while(r ~= M)
    SRBM = zeros(M,N);
    for i=1:N
        ones = randperm(M, S);
        SRBM(ones,i) = 1;
    end
    r = rank(SRBM);
end

for i=1:N
    SRBM(:,i) = SRBM(:,i)*(M/(M+1))^(N-i);
end
    
end

