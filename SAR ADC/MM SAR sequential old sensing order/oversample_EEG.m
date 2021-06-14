dataset = load('/Users/mac/Documents/MATLAB/Thesis/SAR ADC/MM SAR sequential/EEGdata_ch1.mat');
data = dataset.EEGdata_ch1;
output = zeros(3840,80);
for i = 1:80
    data_round = double(data(:,:,i));
    
    range = max(data_round) - min(data_round);
    
    minimum = min(data_round);
    
    data_round = (data_round-minimum)/(range);
   
    data_round = 2*data_round - 1;
    
    data_round = data_round/1000;   %% values are now [-1 mV, 1mV]
    
    data_round = resample(data_round, 3840, 384);
    
    output(:,i) = data_round;
end
save('upsampled_EEG', "output");
