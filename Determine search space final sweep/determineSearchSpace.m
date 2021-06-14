clear;
%% Set constraints
sum=0;
for i=0:10
   sum = sum+2^i; 
end

k = 1.38064852e-23;
T = 293;
Vref = 2;
[MAC_cap, DAC_cap] = getCapSizes(10,3.4878e-09, 1.995e-15, 15, .01, .001025);
minimum_sampling_cap = 12*k*T*2^(2*10)/Vref^2;  
DAC_cap = max([minimum_sampling_cap, 1.995e-15, DAC_cap]);

CUnitMax = sum*ceil(DAC_cap/(1.995e-15)) %max area is area for 14bit ADC
P_max = 3e-6; %SAR power for 6bit with transmission etc


%%  Check search space points
M_array = [75, 150, 192];
N=384;
cap_ratio=15;
number_of_bits_array = 5:7;
leakage_size_boostarray = 1.5:.5:3;
sampling_gain_array = [16, 24, 32];
simulationPoints = [];
simPointAreas = [];
simPointPower = [];
for i = 1:length(M_array)
    for j = 1:length(sampling_gain_array)
        for l = 1:length(number_of_bits_array)
            for x = 1:length(leakage_size_boostarray)
    
    M=M_array(i);
    if (M==75)
        nb_activations = 1;
    else
        nb_activations = 2;
    end
    sampling_gain = sampling_gain_array(j);
    number_of_bits = number_of_bits_array(l);
    leakage_size_boost = leakage_size_boostarray(x);
    
    [power, CUnit] = getSystemCharacteristics(number_of_bits,M, N, nb_activations, cap_ratio, leakage_size_boost, sampling_gain);
    if (power<P_max && CUnit<CUnitMax)
        simulationPoints = [simulationPoints; M, number_of_bits, sampling_gain, leakage_size_boost];
        simPointAreas = [simPointAreas; CUnit];
        simPointPower = [simPointPower; power];
    end
        
    
            end
        end
    end
end

%     fprintf("Estimated sim time: " + num2str(3*size(simulationPoints,1)/24) + " days.\n");
%     fprintf("This is apparently way off. Times vary with chosen M. The 100 points here took about 3 days.\n")
    