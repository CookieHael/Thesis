function ADC_input = database_load(name, epoch)
%DATABASE_LOAD Summary of this function goes here
%   Detailed explanation goes here
EEGdata = load(name);
fields = fieldnames(EEGdata);
data = EEGdata.(fields{1});
selected = data(1,:, epoch);
ADC_input = (2/1000*(selected-min(selected))/(max(selected)-min(selected))) -1/1000;
end

