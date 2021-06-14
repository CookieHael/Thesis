function [capSizeDAC, capSizeMAC] = getCapSizes(NOB,process_param, unitCap, F_per_M)

unitArea = unitCap/F_per_M;

%% DAC
sum=0;
for i=0:NOB-2
    sum = sum + (2^i);
end

%INL
minAreaINL = (6*process_param*sqrt(sum))^2;
capSizeDAC1 = minAreaINL*F_per_M;

%DNL
sum = sum + 2^(NOB-1);
minAreaDNL = (3*process_param*sqrt(sum))^2;
capSizeDAC2 = minAreaDNL*F_per_M;

capSizeDAC = max(capSizeDAC1, capSizeDAC2);

end