function sensingOrder = determineSensingOrder(SRBM, nb_act)
%DETERMINESENSINGORDER Determines in what order channels should be digitzed
%for this specific SRBM

sensingOrder = [];
for i=1:size(SRBM,2)
    c = SRBM(:,i);
    c = (c~=0);
    ind = find(c);
    for j = 1:nb_act
        if (~ismember(ind(j),sensingOrder))
            sensingOrder = [sensingOrder, ind(j)];
        end
    end
end

