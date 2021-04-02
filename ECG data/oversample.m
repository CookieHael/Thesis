directory = dir('/Users/mac/Documents/MATLAB/Thesis/ECG data/MLII/1 NSR/');
for i=3:5
    name=directory(i).name;
    fullname = append('/Users/mac/Documents/MATLAB/Thesis/ECG data/MLII/1 NSR/',name);
    x=load(fullname);
    x=x.val;
    range=max(x)-min(x);
    minim = min(x);
    x=(x-minim)/(range);
    x=2*x - 1;
    x=x/1000;
    x=resample(x, 20000, 360);
    %t=(1/200000)*(1:length(x));
    %y=[t; x];
    %ts = timeseries(y(2:end,:),y(1,:));
    save(convertCharsToStrings(name), "x");
end
