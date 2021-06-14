function ADC_input = database_load(name, epoch)
%DATABASE_LOAD Summary of this function goes here
%   Detailed explanation goes here
EEGdata = load(name);
fields = fieldnames(EEGdata);
data = EEGdata.(fields{1});
selected = squeeze(data);
data = (2/10000*(selected-min(selected))./(max(selected)-min(selected)));
ADC_input = data(:, epoch);

end

