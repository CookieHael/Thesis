function data = load_output_data(str, idx)

f=dir("Output data/"+str+"/*"+".txt");
if str=='C'
    f=dir("Output data/"+str+"/*"+".TXT");
end
name = strcat(f(idx).folder,'/',f(idx).name);
data = readvars(name);
[range, minima] = getRange(str, idx);
data = ((2*((data-minima)/range))/1000)';
end


