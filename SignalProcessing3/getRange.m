function [range, minima] = getRange(input, idx)

f=dir("Neural data/"+input+"_data/*"+".txt");
if input=='C'
    f=dir("Neural data/"+input+"_data/*"+".TXT");
end

name = strcat(f(idx).folder,'/',f(idx).name);
new = readvars(name);
range = max(new)-min(new);
minima = min(new);
end

