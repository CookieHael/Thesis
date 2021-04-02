function [Phi,flag] = genPhi(M,N,D)
% generate the sensing matrix Phi, which is a sparse binary matrix.
%  Input:
%        M - row number of the matrix
%        N - column number of the matrix
%        D - number of nonzero entries of value 1 in each column of Phi
%
% Output
%        Phi - the generated sparse binary matrix
%       flag - if the genrated Phi matrix is full rank, then flag = 1;
%              otherwise, flag= 0.
%
% Zhilin Zhang
%

if M>N,
    error('M should <= N\n');
end

Phi = zeros(M,N);

for i = 1 : N
    ind = randperm(M);
    indice = ind(1:D);
    col = zeros(M,1);
    col(indice) = ones(D,1);
    Phi(:,i) = col;
end

if rank(Phi) == M
    flag = 1;
else
    flag = 0;
end

