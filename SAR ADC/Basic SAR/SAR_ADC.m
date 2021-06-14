clear;

verbose = false;
% for i=1:20
%%  ----------  Circuit parameters  --------------
% number_of_bits = 1+i;
number_of_bits = 6;
input_BW=256;
V_thermal = 25.27*10^-3;
k = 1.38064852e-23;
F_per_M = 0.001025;
T = 293;
sine_amplitude = .001;
sine_bias = 0;
Vref = 2;
noise_quantization = Vref/(sqrt(12)*2^(number_of_bits));
C_mismatch_parameter = 3.4878e-09;


% ADC_input = load('/Users/mac/Documents/MATLAB/Thesis/ECG data/100m (0).mat').x;  %.001 because normalized between -.001 & .001
% ADC_input = ADC_input(2:278);
ADC_input = database_load('EEGdata_ch1.mat', 1)'/10;

% Clock & sine timing settings
% sine_frequency_rad = 2*pi*input_BW;  % Set to 256Hz
%nyquist_frequency = sine_frequency*number_of_bits*2;
clk_freq = 2.1*(number_of_bits+1)*input_BW;
clk_period = 1/clk_freq;

%% Circuit non-idealities

non_idealities_on = 1;
% Sample time jitter
sampling_jitter_on = non_idealities_on;
sample_time_jitter = 16e-12;

% Cap sizing
minimum_technology_cap = 1.995e-15;
minimum_sampling_cap = 12*k*T*2/3*2^(2*number_of_bits)/Vref^2;
minimum_matching_cap = getCapSizes(number_of_bits,C_mismatch_parameter, minimum_technology_cap, F_per_M);
minimum_cap = max([minimum_technology_cap,minimum_sampling_cap, minimum_matching_cap]);

% Switching noise
switch_noise_on = non_idealities_on;
Cs = minimum_cap*10;
Cf = minimum_cap;
sampling_gain=Cs/Cf;
switching_noise = k*T/Cf;


% Comparator offset/hysteresis & noise
comparator_offset_on = non_idealities_on;
comparator_high_offset = 0.2*10e-4;
comparator_low_offset = -.1*10e-4;
comparator_noise_on = non_idealities_on;
comparator_noise_rms = 1*k*T* 2/3 / minimum_cap;



%% LNA
LNA_noise_on = non_idealities_on;
LNA_gain = 1000;
LNA_bandwidth = 5*input_BW; %% in Hz
LNA_NEF = 1.08;
gmoverid = 20; %%I=150nA

% Determine minimum current if bandwidth-constrained
gm_min_1 = LNA_bandwidth*2*pi*minimum_cap;
Id_min_1 = gm_min_1*gmoverid;
LNA_input_referred_noise_1 = LNA_NEF/(sqrt(2*Id_min_1/(pi*V_thermal*4*k*T*LNA_bandwidth)));

% Determine minimum current if slewrate-constrained
SR_required = Vref*clk_freq; %Minimum current
Id_min_2 = SR_required*minimum_cap;
gm_min_2 = Id_min_2*gmoverid;
LNA_input_referred_noise_2 = LNA_NEF/(sqrt(2*Id_min_2/(pi*V_thermal*4*k*T*LNA_bandwidth)));

% Determine minimum current if noise-constrained
LNA_input_referred_noise_3 = noise_quantization/(LNA_gain*sampling_gain);
Id_min_3 = (LNA_NEF/LNA_input_referred_noise_3)^2*pi*4*k*T*LNA_bandwidth*V_thermal;

% Determine which of the three cases limits LNA performance
I_LNA = max([Id_min_1, Id_min_2, Id_min_3]);
if (I_LNA==Id_min_3 && Id_min_3>130e-09)
    warning("gm/Id not really valid anymore, 20 at 130nA but lower at higher currents");
end

LNA_input_referred_noise = min([LNA_input_referred_noise_1,LNA_input_referred_noise_2,LNA_input_referred_noise_3]); % Take maximum current of 3 required, noise will be minimum due to NEF formula
LNA_SR = I_LNA/(minimum_cap);

if(sine_amplitude>Vref/(2*sampling_gain*LNA_gain))
    warning("Careful! Input signal might be saturating ADC conversion range.\n");
end

%% Compression    no longer implemented

% % ZOH compression
% ZOH_delta = 100*Vref/(2^number_of_bits);
% 
% % Truncation compression
% 
% number_of_bits_transmission = 4;

%%% ----------------  Testing  ----------------

%% Simulate
sim_out = sim('ADC', length(ADC_input)*(number_of_bits+1)*clk_period);


%% SNDR calculation
% sndr = sinad_ADC(sim_out.yout.signals(1).values(3:end), number_of_bits, clk_freq);
% rms_out = calculate_rms(1000*reshape_output(squeeze(sim_out.yout.signals(2).values), number_of_bits), reshape_output(sim_out.yout.signals(1).values(2:end), number_of_bits));
% rms_out = calculate_rms(1000*ADC_input(2:end-1), sim_out.yout.signals(1).values(3:end));

% rms_out
% sndr_minus_ideal = sndr - 6.02*number_of_bits - 1.76 + 20*log(Vref/((2*sine_amplitude)*LNA_gain*sampling_gain));

%% Power calculations
C_0_DAC = minimum_cap;
P_dac = (2^number_of_bits*clk_freq*C_0_DAC)/(number_of_bits+1) * ( ((5/6) - (1/2)^number_of_bits -(1/3)*(1/2)^(2*number_of_bits))*Vref^2 - 1/2*(sine_bias)^2 - (1/2)^number_of_bits*Vref*(sine_bias));

V_eff = 0.1;
V_fs = Vref;
C_load_sar = minimum_cap;
P_comp = 2*number_of_bits*log(2)*clk_freq*number_of_bits/(number_of_bits+1)*V_eff*V_fs*C_load_sar;

alpha_logic = 0.4;
C_logic = .77*1e-15;
P_logic = alpha_logic * (2*number_of_bits+1)*8*C_logic*Vref^2*(clk_freq-clk_freq/(number_of_bits+1));

P_tot_SAR = (P_dac + P_comp + P_logic);

V_thermal = 25.27*10^-3;
I_LNA = (LNA_NEF/LNA_input_referred_noise)^2*pi*4*k*T*LNA_bandwidth*V_thermal;
P_LNA = Vref*I_LNA;

I_SH = 2*sampling_gain*C_logic*Vref*(clk_freq/(number_of_bits+1));
P_SH = Vref * I_SH;

% Transmission power

transmission_power_per_bit = 10^-9;
P_transmission_full_signal = clk_freq/(number_of_bits+1) * number_of_bits * transmission_power_per_bit;


P_tot = P_tot_SAR+P_LNA + P_transmission_full_signal + P_SH;
P_LNA_percentage = 100*P_LNA/ (P_tot);
P_SAR_percentage =  100*P_tot_SAR/ (P_tot);
P_comp_percentage = 100*P_comp/P_tot;
P_dac_percentage = 100*P_dac/P_tot;
P_SH_percentage = 100*P_SH/P_tot;
P_logic_percentage = 100*P_logic/P_tot;
P_transmission_full_percentage = 100*P_transmission_full_signal/P_tot;

figure
hold on
plot(10000*ADC_input(1:end-1))
plot(sim_out.yout.signals(1).values(3:end))
hold off
xlabel("Sample index")
ylabel("Output voltage (V)")
legend("Input signal", "Output signal")

20*log10(10000*ADC_input(1:end-1)/(10000*ADC_input(1:end-1)-sim_out.yout.signals(1).values(3:end)'))


% figure
% hold on
% bar([P_LNA_percentage, P_SH_percentage; P_comp_percentage; P_dac_percentage;P_logic_percentage]);
% set(gca, 'XTick', 1:5, 'XTickLabels', { "LNA", "S&H", "Comparator", "DAC", "Logic"});
% ylabel("Relative power draw percentages without tranmission")
% title("Total power: " + P_tot*1e6 + "\muW")
% axes('Position',[.7 .7 .2 .2])
% box on
% bar([P_LNA_percentage; P_SAR_percentage; P_transmission_full_percentage; ]);
% set(gca, 'XTick', 1:3, 'XTickLabels', {"LNA", "SAR", "Transmission"});
% title("With TX")

figure
hold on
bar([P_SH_percentage; P_comp_percentage; P_dac_percentage;P_logic_percentage]);
set(gca, 'XTick', 1:4, 'XTickLabels', { "S&H", "Comparator", "DAC", "Logic"});
ylabel("Relative power draw percentages without Tx & LNA")
title("Total power: " + P_tot*1e6 + "\muW")
axes('Position',[.7 .7 .2 .2])
box on
bar([P_LNA_percentage; P_SAR_percentage; P_transmission_full_percentage; ]);
set(gca, 'XTick', 1:3, 'XTickLabels', {"LNA", "SAR", "Transmission"});
title("With TX and LNA")



% Compressed signals
% 
%  sndr_compression1 = sinad_ADC(sim_out.yout.signals(2).values, number_of_bits, clk_freq);
%  P_transmission_compression1 = clk_freq/(number_of_bits+1) * number_of_bits_transmission * transmission_power_per_bit;
%  compression1_ratio = round(100*(number_of_bits-number_of_bits_transmission)/number_of_bits,2);   %Voorstelling? Relatief vs absoluut
%  P_tot_compression1 = P_LNA+P_tot_SAR + P_transmission_compression1;

% sndr_ZOHcompression = sinad_ADC(sim_out.yout.signals(4).values, number_of_bits, clk_freq);
% P_transmission_ZOH = (1/(clk_period*20000/sim_out.yout.signals(3).values(size(sim_out.yout.signals(3).values,1))))* number_of_bits * transmission_power_per_bit;


%% Print results numbers
% fprintf("ADC resolution = " + number_of_bits + " bits.\n");
% fprintf("Calculated SNDR is " + sndr + " dB. This is " + sndr_minus_ideal + "dB from ideal. \n");
% if (verbose)
%     fprintf("Quantization noise = " + 10^6 * noise_quantization + "uVrms.\n");
%     fprintf("Minimal capacitor value = " + round(minimum_sampling_cap,4) + "fF.\n");
%     fprintf("Chosen sampling capacitor =  "+ Cs*10^9 +"nF. If larger than minimum, this lowers ADC power efficiency.\n")
% end
% fprintf("Total power draw =  " + P_tot + " nW.\n");
% if (verbose)
%     fprintf("LNA power draw = " + P_LNA + " nW, or " + P_LNA_percentage + "%% of the total power draw.\n");
%     fprintf("Comparator power draw = " + P_comp + " nW, or " + P_comp_percentage + "%% of the total power draw.\n");
%     fprintf("DAC power draw = " + P_dac + " nW, or " + P_dac_percentage + "%% of the total power draw.\n");
% end
% fprintf("Transmission power draw = " + P_transmission_full_signal + " nW, or " + P_transmission_full_percentage + "%% of the total power draw.\n");
% fprintf("Truncated signal's SNDR = " + sndr_compression1 + " dB, with compression ratio of " + compression1_ratio + "%%.\nTotal compression 1 power draw = " + P_tot_compression1 + "nW, or " + round(100*P_tot_compression1/P_tot, 1) + "%% of uncompressed total power.\n");
% end

%% Sample time jitter sweep

% sampling_jitter_on = 1;
% sample_time_jitters = [0 linspace(1e-5, 1e-8, 50)];
% 
% SNDR = zeros(size(sample_time_jitters));
% 
% for j = 1:size(sample_time_jitters,2)
%     sample_time_jitter = sample_time_jitters(j);
%     sim_out = sim('ADC', 20000*clk_period);
%     SNDR(j) = sinad_ADC(sim_out.yout.signals.values, number_of_bits, clk_freq);
% end
% f=gcf;
% scatter(sample_time_jitters, SNDR)
% exportgraphics(f, 'sampling_time_jitter_sweep.pdf');
% sampling_jitter_on = 0;


%% Switching noise sweep

% switch_noise_on = 1;
% k = 1.38064852e-23;
% T = 293;
% nb_points = 50;
% Cs_array = linspace(1e-15, 1e-14, nb_points);
% Cf_array = Cs_array/sampling_gain;
% nb_points = size(Cs_array, 2);
% SNDR = zeros(1,nb_points);
% 
% for j = 1:nb_points
%     Cs = Cs_array(j);
%     Cf = Cf_array(j);
%     sampling_gain=Cs/Cf;
%     sim_out = sim('ADC', 20000*clk_period);
%     SNDR(j) = sinad_ADC(sim_out.yout.signals.values, number_of_bits, clk_freq);
% end
% f=gcf;
% scatter(Cs_array, SNDR);
% exportgraphics(f, 'switching_noise_sweep.pdf');
% switch_noise_on = 0;


%% Comparator offset

% comparator_offset_on = 1;
% comparator_high_offset_array = linspace(1e-5, 5e-1, 50);
% comparator_low_offset_array = comparator_high_offset_array/2;
% nb_points = size(comparator_high_offset_array,2);
% SNDR = zeros(1,nb_points);
% 
% for j = 1:nb_points
%     comparator_high_offset = comparator_high_offset_array(j);
%     comparator_low_offset = comparator_low_offset_array(j);
%     sim_out = sim('ADC', 20000*clk_period);
%     SNDR(j) = sinad_ADC(sim_out.yout.signals.values, number_of_bits, clk_freq);
% end
% f=gcf;
% scatter(comparator_high_offset_array, SNDR);
% exportgraphics(f, 'comparator_offset_sweep.pdf');
% comparator_offset_on = 0;

% %% LNA noise
% 
% nb_points = 50;
% LNA_noise_on = 1;
% LNA_bandwidth = 10000;
% LNA_array = linspace(0, 5*1e-3, nb_points);
% SNDR = zeros(1,nb_points);
% 
% for j = 1:nb_points
%     LNA_input_reffered_noise = LNA_array(j);
%     sim_out = sim('ADC', 20000*clk_period);
%     SNDR(j) = sinad_ADC(sim_out.yout.signals.values, number_of_bits, clk_freq);
% end
% f=gcf;
% scatter(comparator_high_offset_array, SNDR);
% exportgraphics(f, 'LNA_noise_sweep.pdf');
% LNA_noise_on = 0;
%     

2^number_of_bits*C_0_DAC/(1.995e-15)

%% Helper functions
function sndr_out = sinad_ADC(input_vector, number_of_bits, clk_freq)
% y = zeros(round((size(input_vector,1))/(2*(number_of_bits+1))),1);
% i=1;
% for j = 1:size(y,1)
%     y(j) = input_vector(i);
%     i= i + 2 + 2*number_of_bits;
% end
%sinad(y, clk_freq/(number_of_bits+1));
sndr_out = sinad(input_vector, clk_freq/(number_of_bits+1));
end

function out = reshape_output(input, number_of_bits)
out = zeros(round((size(input,1))/(2*(number_of_bits+1))),1);
i=1;
    for j = 1:size(out,1)
        out(j) = input(i);
        i= i + 2 + 2*number_of_bits;
    end
end

function rms_out = calculate_rms(input_values, output_values)
hold on
plot(input_values);
plot(output_values);
legend("Input", "Output");
hold off

if (length(input_values) ~= length(output_values))
    rms_out = -inf;
    return
end

rms_out = sqrt(sum( (input_values-output_values).^2 )  / length(output_values) );
end

