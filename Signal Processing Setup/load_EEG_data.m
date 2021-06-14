function data = load_EEG_data(str, idx)

f=dir("Neural data/"+str+"_data/*"+".txt");
if str=='C'
    f=dir("Neural data/"+str+"_data/*"+".TXT");
end
name = strcat(f(idx).folder,'/',f(idx).name);
data = readvars(name);
data = data(1:4096);
data = ((2*((data-min(data))/(max(data)-min(data))))/1000)';
end

