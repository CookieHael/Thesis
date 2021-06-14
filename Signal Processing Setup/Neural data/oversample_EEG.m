    data_round = data(1:floor(length(data)/10), i);
    
    range = max(data_round) - min(data_round);
    
    minimum = min(data_round);
    
    data_round = (data_round-minimum)/(range);
   
    data_round = 2*data_round - 1;
    
    data_round = data_round/1000;   %% values are now [-1 mV, 1mV]
    
    data_round = resample(data_round, 8000, 256);
    
    output(:,i) = data_round;

