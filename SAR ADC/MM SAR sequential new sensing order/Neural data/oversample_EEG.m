dataset = load('/Users/mac/Documents/MATLAB/Thesis/SAR ADC/MM SAR/Neural data/A01S.mat');
data = dataset.data{1,1}.X;
output = zeros(156813, 16);
for i = 1:16
    data_round = data(1:floor(length(data)/10), i);
    
    range = max(data_round) - min(data_round);
    
    minimum = min(data_round);
    
    data_round = (data_round-minimum)/(range);
   
    data_round = 2*data_round - 1;
    
    data_round = data_round/1000;   %% values are now [-1 mV, 1mV]
    
    data_round = resample(data_round, 8000, 256);
    
    output(:,i) = data_round;
end
save('upsampled_EEG', "output");
