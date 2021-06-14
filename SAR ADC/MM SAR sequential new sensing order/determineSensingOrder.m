function [sensingMatrix, sensingOrder] = determineSensingOrder(SRBM, nb_act, nb_ADC)
%DETERMINESENSINGORDER Determines in what order channels should be digitzed
%for this specific SRBM. Also, create a matrix with 0's until the cycle in
%which a specific cap is digitized and ones afterwards.

% sensingOrder = [];
% 
% 
% lastAcc = zeros(1,size(SRBM,1));
% for i = 1:size(SRBM,1)   % Find last accumulation on each row
%     row = SRBM(i, :);
%     ind = find(row);
%     lastAcc(1,i) = ind(end);
% end
% 
% lastAcc = lastAcc-min(lastAcc);   % Normalize wrt first cycle smth should be digitized

sensingOrder_prep = [];
for i=size(SRBM,2):-1:1
    c = SRBM(:,i);
    c = (c~=0);
    ind = find(c);
    for j = 1:nb_act
        if (~ismember(ind(j),sensingOrder_prep))
            sensingOrder_prep = [sensingOrder_prep, ind(j)];
        end
    end
end

if size(sensingOrder_prep)<size(SRBM, 1)
    for i=1:size(SRBM, 1)
        if ~ismember(i,sensingOrder_prep)
             sensingOrder_prep = [sensingOrder_prep, i];
        end
    end
end

sensingOrder = [];
sensingOrder_prep = flip(sensingOrder_prep); %Calculation above gives reverse sensing orders

for i=1:nb_ADC
    sensingOrder = [sensingOrder,sensingOrder_prep(:,i:nb_ADC:end)];
end

sensingMatrix = zeros(size(SRBM));
for i = 1:size(SRBM,1)
    row = SRBM(i,:);
    idx = find(row);
    if (~isempty(idx))
        if (idx(end)<size(SRBM,2))
            sensingMatrix(i, idx(end)+1:size(SRBM,2)) = 1;
        else
            sensingMatrix(i, idx(end):size(SRBM,2)) = 1;
        end
    end
end
