function sensingOrder = determineSensingOrder(SRBM, nb_act, nb_ADC)
%DETERMINESENSINGORDER Determines in what order channels should be digitzed
%for this specific SRBM

sensingOrder_prep = [];
for i=1:size(SRBM,2)
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

for i=1:nb_ADC
    sensingOrder = [sensingOrder,sensingOrder_prep(:,i:nb_ADC:end)];
end

