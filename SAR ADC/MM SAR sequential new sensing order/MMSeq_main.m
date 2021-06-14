clear;
verbose = false;
% tic

%%  Constants &  process parameters
k = 1.38064852e-23;
T = 293;
V_thermal = 25.27*10^-3;
minimum_technology_cap = 1.995 * 1e-15;
F_per_M = 0.001025;
area_per_C = 1/ F_per_M;  %1/ (F/m^2)
C_mismatch_parameter = 3.4878e-09;
C_logic = 1*1e-15;
noise_multiplier = 1;
Vref = 2;
gmoverid = 20; %%I=150nA
V_eff = 1/gmoverid;
fmax_ADC = 10e6;   %% Determine from simulation!


%% Inputs
epoch = 1;
rsnr_uncorr = zeros(80,1);
ssim_uncorr = zeros(80,1);
nmse_0 = zeros(80,1);
ssim_0 = zeros(80,1);
snr_0 = zeros(80,1);
nmse_1 = zeros(80,1);
ssim_1 = zeros(80,1);
snr_1 = zeros(80,1);
nmse_2 = zeros(80,1);
ssim_2 = zeros(80,1);
snr_2 = zeros(80,1);
nmse_smooth = zeros(80,1);
ssim_smooth = zeros(80,1);
snr_smooth = zeros(80,1);
tempRSNR = 0;
for p=1:1
    p
ADC_input = database_load('EEGdata_ch1.mat',p)';

% ADC_input_OS =  database_load('upsampled_EEG', epoch)'/10;
%data = load('100m (0).mat').x/10;


%% System settings
number_of_bits = 6;
input_BW = 256; %This should be the Nyquist frequency, and input should be oversampled wrt this
signal_peak_amplitude = .00002;
signal_bias = 0;
nb_MAC = 384;
cap_ratio = 15;
nb_channels = 150;
tertiary_matrix = false;
nb_activations = 1;
cap_gain_error = 1e-2;

non_idealities_on = 1;
noise_on = 1;
leakage_on = 1;

%% Warnings
if ~noise_on
    warning("Noise currently off!")
end

if noise_multiplier>1 
    warning("Noise multiplier not one!")
end

if ~leakage_on
    warning("Leakage currently  off!")
end

noise_quantization = sqrt((Vref)^2/12*2^(-2*number_of_bits));

%% Sensing matrix setup

if(nb_channels == 150 && nb_MAC == 384 && nb_activations == 2)
    sensing_matrix_large = load("sensing150.mat").sensing_matrix_large;  %%Capacitor ratio was 35
    sensing_matrix_corrected = load("sensing150.mat").sensing_matrix_corrected;
    
elseif (nb_channels == 192 && nb_MAC == 384 && nb_activations == 2)
    sensing_matrix_large = load("sensing192.mat").Phi;
    sensing_matrix_corrected = load("sensing192.mat").corrected_Phi;
elseif (nb_channels == 50 && nb_MAC == 384 && nb_activations == 1)
    sensing_matrix_large = load("Phi50ch384MAC.mat").phi;
    sensing_matrix_corrected =  rescaleSRBM(sensing_matrix_large, cap_ratio);
elseif (nb_channels == 75 && nb_MAC == 384 && nb_activations == 1)
    sensing_matrix_large = load("phi384MAC75ch.mat").bestPhi;
    sensing_matrix_corrected =  rescaleSRBM(sensing_matrix_large, cap_ratio);
else
    [sensing_matrix_large, sensing_matrix_corrected] = generateCorrectedSRBM(nb_channels, nb_MAC, nb_activations, cap_ratio);
end
if (cap_ratio ~=35)
    sensing_matrix_corrected = rescaleSRBM(sensing_matrix_large, cap_ratio);
end

if( nb_MAC ~= size(sensing_matrix_large, 2) || nb_channels ~= size(sensing_matrix_large,1))
    error("Wrong matrix dimensions");
end
if (tertiary_matrix)
    [sensing_matrix_corrected, sensing_matrix_large] = makeTertiary(sensing_matrix_corrected, sensing_matrix_large);
end

%% Clock & sine timing settings

clk_freq = 2.1*input_BW;
clk_period = 1/clk_freq;

%% Sample time jitter
sampling_jitter_on = noise_on;
sample_time_jitter = 16e-12;

%% Capacitor sizing
minimum_sampling_cap = 2 * 12*k*T*2^(2*number_of_bits)/Vref^2;  % Double to increase noise performance
[minimum_match_cap_DAC, minimum_match_cap_MAC] = getCapSizes(number_of_bits,C_mismatch_parameter, minimum_technology_cap, cap_ratio, cap_gain_error, F_per_M);
MAC_cap = max([minimum_sampling_cap, minimum_technology_cap, minimum_match_cap_MAC]);
DAC_cap = max([minimum_sampling_cap, minimum_technology_cap, minimum_match_cap_DAC]);

%% MAC non-idealities & noise
C_sense_size_multiplier = 2.5;
C1 = C_sense_size_multiplier*MAC_cap;
C2 = C1*cap_ratio; 
MM_incomplete_transfer_coeff = C1/(C1+C2);
if (non_idealities_on)
    MM_attenuation_coeff = C2/(C1+C2);  %Badly chosen name: this is voltage attenuation by charge sharing architecture
else
    MM_attenuation_coeff = 1;  %Actually more like 1-leakage
end
MM_sampling_noise = noise_multiplier*noise_on*((k*T * C1/(C1+C2)^2) + (k*T*(C1*C2/(C1+C2))/C2^2));  % Formula with square root is sigma, model uses variance

%% Leakage
leak_current = leakage_on *1*10^(-15);

%% Sample & hold (with gain)
switch_noise_on = noise_on;
Cs_unit_cap = 4/(MM_incomplete_transfer_coeff);
Cf_unit_cap = 2;
sampling_gain = Cs_unit_cap/Cf_unit_cap;
switching_noise = noise_multiplier*k*T/(Cf_unit_cap*minimum_sampling_cap);


%% Comparator offset/hysteresis & noise
comparator_offset_on = 1;
comparator_high_offset = 655e-6;
comparator_low_offset = comparator_high_offset;
comparator_noise_on = 1;
comparator_noise_rms = noise_multiplier* 1*sqrt(k*T/C_logic * (2/3)) ;


%% LNA
LNA_noise_on = 1;
LNA_gain = 1000;
LNA_bandwidth = 3*input_BW; %% in Hz
LNA_NEF = 1.08;

% Determine minimum current if bandwidth-constrained
gm_min_1 = LNA_gain*LNA_bandwidth*2*pi*nb_activations*MAC_cap;
Id_min_1 = gm_min_1*gmoverid;
LNA_input_referred_noise_1 = LNA_NEF/(sqrt(2*Id_min_1/(pi*V_thermal*4*k*T*LNA_bandwidth)));

% Determine minimum current if slewrate-constrained
SR_required = Vref*clk_freq; %Minimum current
Id_min_2 = SR_required*nb_activations*MAC_cap;
LNA_input_referred_noise_2 = LNA_NEF/(sqrt(2*Id_min_2/(pi*V_thermal*4*k*T*LNA_bandwidth)));

% Determine minimum current if noise-constrained
LNA_input_referred_noise_3 = 0.5*noise_quantization/(LNA_gain*sampling_gain*MM_incomplete_transfer_coeff);
Id_min_3 = (LNA_NEF/LNA_input_referred_noise_3)^2*pi*4*k*T*LNA_bandwidth*V_thermal;

% Determine which of the three cases limits LNA performance
I_LNA = max([Id_min_1, Id_min_2, Id_min_3]);
if (I_LNA==Id_min_1 && Id_min_3>130e-09)
    warning("gm/Id not really valid anymore, 20 at 130nA but lower at higher currents");
end

LNA_in_rms = min([LNA_input_referred_noise_1,LNA_input_referred_noise_2,LNA_input_referred_noise_3]); % Take maximum current of 3 required, noise will be minimum due to NEF formula
LNA_SR = I_LNA/(nb_activations*MAC_cap);
pd = makedist('Normal', 0, LNA_in_rms);
rng(420699);
LNA_noise_vector = [pd.random(nb_MAC,1);0]';


%% Total gain check
total_gain = LNA_gain*sampling_gain*MM_incomplete_transfer_coeff;
if(signal_peak_amplitude*total_gain>Vref/2)
    warning("Careful! Input signal might be saturating ADC conversion range.\n");
end


%% ZOH settings

ZOH_delta = 0*Vref/(2^number_of_bits);
if (ZOH_delta>0)
    fprintf("Reminder: ZOH is on, signal quality may be degraded.\n");
end


%% Determine minimum required number of ADCs

OCR = nb_activations*(number_of_bits+3);
nb_ADC = ceil(OCR*clk_freq/fmax_ADC);


%% Simulate

[sensingOrderMatrixFull, sensing_order] = determineSensingOrder(sensing_matrix_large, nb_activations, nb_ADC);
sensingOrderMatrix = sensingOrderMatrixFull(sensing_order,:);
sensingOrderMatrix = sensingOrderMatrix'; %Transpose because next operation puts columns under eachother, we want rows appended instead
sensingOrderMatrix = sensingOrderMatrix(:);
sensing_matrix_sensing_order = sensing_matrix_large(sensing_order,:);
sensing_matrix_prep = sensing_matrix_sensing_order';
sensing_matrix = sensing_matrix_prep(:)';


sim_number = randi(1e6); %Set this to fixed number for repeatable results. In thesis, 10 was used
sim_noise_seed = sim_number;


sim_out = sim('MMSAR_Sequential', (nb_channels+1)*clk_period*(nb_MAC+1)+clk_period*(number_of_bits+5));

cs_out = zeros(nb_channels,1);
for i = 1:nb_channels
    cs_out(sensing_order(i)) = sim_out.yout.signals(1).values(i+3);
end
nb_ZOH = 0*sim_out.yout.signals(2).values(end);


%% Compressive sensing reconstruction      
if (non_idealities_on)
    A = sensing_matrix_corrected*wmpdictionary(nb_MAC, 'lstcpt', {'dct'});
else
    A = sensing_matrix_large*wmpdictionary(nb_MAC, 'lstcpt', {'dct'});
end

recovery = BSBL_BO(A, cs_out, 1:13:nb_MAC, 2, 'prune_gamma',-1, 'max_iters',15);
recovered_signal = (idct(recovery.x))';
recovered_signal = recovered_signal(1:nb_MAC);

rsnr_uncorr(p) = getRSNR(ADC_input*total_gain, recovered_signal);
windowLen = 100;
[mssim, ssim_map] = ssim_1d(recovered_signal, ADC_input*total_gain, windowLen);
ssim_uncorr(p) = mssim;

%   RMS calculations
expected_coeff = sensing_matrix_corrected*(total_gain*ADC_input');
expected_recovery = BSBL_BO(A, expected_coeff, 1:13:nb_MAC, 2, 'prune_gamma',-1, 'max_iters',15);
expected_recovered_signal = (idct(expected_recovery.x))';

[rms_out, mse_out] = calculateRMS(ADC_input*total_gain, recovered_signal);
RSNR_out_expected = getRSNR(ADC_input*total_gain, expected_recovered_signal);
[rms_out_expected, mse_out_expected] = calculateRMS(expected_recovered_signal, ADC_input*total_gain);

% SSIM calculations
windowLen = 100;
[mssim, ssim_map] = ssim_1d(recovered_signal, ADC_input*total_gain, windowLen);
% ssim_1 = mssim;

[mssim_expected, ssim_map_expected] = ssim_1d(ADC_input*total_gain, expected_recovered_signal, windowLen);
ssim_expected = mssim_expected;


%% Correct leakage
cs_out_LC = negateLeakage(sensing_matrix_large, sensingOrderMatrixFull, sampling_gain, leak_current, C2, clk_period, cs_out);
recovery = BSBL_BO(A, cs_out_LC, 1:13:nb_MAC, 2, 'prune_gamma',-1, 'max_iters',15);

recovered_signal_LC = (idct(recovery.x))';
recovered_signal_LC = recovered_signal_LC(1:nb_MAC);
[rms_out_LC, mse_out_LC] = calculateRMS(ADC_input*total_gain, recovered_signal_LC);
mse_out_LC;
RSNR_out_LC = getRSNR(ADC_input*total_gain, recovered_signal_LC);
RSNR_out_LC

[mssim, ssim_map] = ssim_1d(recovered_signal_LC, ADC_input*total_gain, windowLen);
ssim_1_LC = mssim
if ssim_1_LC>tempRSNR
    temp=sensing_matrix_large;
    tempRSNR=ssim_1_LC;
end
% %% Second reconstruction method (set for sparse signals, this sometimes performs better for high leakage&compression, test out) 
recovery = BSBL_BO(A, cs_out_LC, 1:13:nb_MAC, 2);
recovered_signal_2 = (idct(recovery.x))';
recovered_signal_2 = recovered_signal_2(1:nb_MAC);
[rms_out_2, mse_out_2] = calculateRMS(ADC_input*total_gain, recovered_signal_2);
nmse_smooth(p) = mse_out_2;
[mssim, ssim_map] = ssim_1d(recovered_signal_2, ADC_input*total_gain, windowLen);
ssim_smooth(p) = mssim;
snr_smooth(p) = getRSNR(ADC_input*total_gain, recovered_signal_2);


%%Try setting 0
cs_out_LC = negateLeakage(sensing_matrix_large, sensingOrderMatrixFull, sampling_gain, leak_current, C2, clk_period, cs_out);
recovery = BSBL_BO(A, cs_out_LC, 1:13:nb_MAC, 0, 'prune_gamma',-1, 'max_iters',15);
recovered_signal_LC = (idct(recovery.x))';
recovered_signal_LC = recovered_signal_LC(1:nb_MAC);
[rms_out_LC, mse_out_LC] = calculateRMS(ADC_input*total_gain, recovered_signal_LC);
mse_out_LC;
RSNR_out_LC = getRSNR(ADC_input*total_gain, recovered_signal_LC);
RSNR_out_LC;
[mssim, ssim_map] = ssim_1d(recovered_signal_LC, ADC_input*total_gain, windowLen);

nmse_0(p) = mse_out_LC;
ssim_0(p) = mssim;
snr_0(p) = RSNR_out_LC;

%%Try setting 1
cs_out_LC = negateLeakage(sensing_matrix_large, sensingOrderMatrixFull, sampling_gain, leak_current, C2, clk_period, cs_out);
recovery = BSBL_BO(A, cs_out_LC, 1:13:nb_MAC, 1, 'prune_gamma',-1, 'max_iters',15);
recovered_signal_LC = (idct(recovery.x))';
recovered_signal_LC = recovered_signal_LC(1:nb_MAC);
[rms_out_LC, mse_out_LC] = calculateRMS(ADC_input*total_gain, recovered_signal_LC);
mse_out_LC;
RSNR_out_LC = getRSNR(ADC_input*total_gain, recovered_signal_LC);
RSNR_out_LC;
[mssim, ssim_map] = ssim_1d(recovered_signal_LC, ADC_input*total_gain, windowLen);

nmse_1(p) = mse_out_LC;
ssim_1(p) = mssim;
snr_1(p) = RSNR_out_LC;
%%Try out of some other setting
cs_out_LC = negateLeakage(sensing_matrix_large, sensingOrderMatrixFull, sampling_gain, leak_current, C2, clk_period, cs_out);
recovery = BSBL_BO(A, cs_out_LC, 1:13:nb_MAC, 2, 'prune_gamma',-1, 'max_iters',15);
recovered_signal_LC = (idct(recovery.x))';
recovered_signal_LC = recovered_signal_LC(1:nb_MAC);
[rms_out_LC, mse_out_LC] = calculateRMS(ADC_input*total_gain, recovered_signal_LC);
RSNR_out_LC = getRSNR(ADC_input*total_gain, recovered_signal_LC);
[mssim, ssim_map] = ssim_1d(recovered_signal_LC, ADC_input*total_gain, windowLen);


nmse_2(p) = mse_out_LC;
ssim_2(p) = mssim;
snr_2(p) = RSNR_out_LC;
end

%% Power calculations with ZOH
% C_0_DAC = DAC_cap;
% E_dac = (2^(number_of_bits) * C_0_DAC)/(number_of_bits+3) * ( ((5/6) - (1/2)^number_of_bits -(1/3)*(1/2)^(2*number_of_bits))*Vref^2 - 1/2*(signal_bias)^2 - (1/2)^number_of_bits*Vref*(signal_bias));
% P_dac = E_dac * nb_channels * nb_activations*(number_of_bits+3) * clk_freq/(nb_MAC) * (nb_channels-nb_ZOH);
% 
% alpha_logic = 0.4;
% P_logic_SAR = alpha_logic * ( (2*number_of_bits+5)*8*C_logic) * Vref^2 * nb_channels/(nb_MAC)*clk_freq*number_of_bits/(number_of_bits+3) * nb_activations*(number_of_bits+3); %include leakage power  2n+1=> 2n+5 (2 extra FF for additional cycles, 2 extra for delaying outputs), + 3C_logic for ZOH gates
% P_logic_ZOH = alpha_logic * C_logic * 17 * Vref^2 * nb_channels/(nb_MAC)*clk_freq * 2/(number_of_bits+3) * nb_activations*(number_of_bits+3); % 3 * AND + 1 * OR + 1 * NOT assuming AND&OR are 4*C_logic
% P_logic_phi = (log2(nb_MAC) + 1) * nb_MAC * 8*C_logic * Vref^2 * clk_freq;  %alpha logic assumed 1 for shift register
% P_logic = P_logic_SAR + P_logic_ZOH + P_logic_phi;
% 
% P_switch = 4*C_logic* (nb_MAC+nb_activations)* Vref^2*clk_freq;
% 
% V_fs = Vref;
% C_load_sar = 4*C_logic;  %Assume latch used for input to SAR
% E_comp = 2*number_of_bits*log(2)*V_eff*V_fs*C_load_sar;
% P_comp = nb_channels/(nb_MAC) * clk_freq * (number_of_bits+3) * ( (nb_channels-nb_ZOH)/nb_MAC * E_comp + 3/(number_of_bits+3) * (1-(nb_channels-nb_ZOH)/nb_MAC) * E_comp);
% 
% P_tot_SAR = (P_dac+P_comp+P_logic);
% 
% P_LNA = Vref*I_LNA;
% 
% I_SH = 2*sampling_gain*minimum_sampling_cap*Vref*clk_freq*nb_channels/(nb_MAC);
% P_SH = Vref * I_SH;
% 
% Transmission power
% 
% transmission_power_per_bit = 10^-9;
% P_transmission_full_signal = clk_freq/nb_MAC * (nb_channels-nb_ZOH) * number_of_bits * transmission_power_per_bit;
% 
% 
% P_tot = P_tot_SAR + P_LNA + P_SH + P_transmission_full_signal+P_switch;
% P_LNA_percentage = 100*P_LNA/P_tot;
% P_SH_percentage = 100*P_SH/P_tot;
% P_comp_percentage = 100*P_comp/P_tot;
% P_dac_percentage = 100*P_dac/P_tot;
% P_SAR_percentage = 100*P_tot_SAR/P_tot;
% P_switch_percentage = 100*P_switch/P_tot;
% P_logic_percentage = 100*P_logic/P_tot;
% P_transmission_full_percentage = 100*P_transmission_full_signal/P_tot;

%% Power calculations without ZOH
C_0_DAC = DAC_cap;
E_dac = (2^(number_of_bits) * C_0_DAC)/(number_of_bits+3) * ( ((5/6) - (1/2)^number_of_bits -(1/3)*(1/2)^(2*number_of_bits))*Vref^2 - 1/2*(signal_bias)^2 - (1/2)^number_of_bits*Vref*(signal_bias));
P_dac = E_dac * nb_channels*clk_freq/(nb_MAC) * (nb_channels);

alpha_logic = 0.4;
P_logic_SAR = alpha_logic * (2*number_of_bits*8*C_logic) * Vref^2 * nb_channels/(nb_MAC)*clk_freq*number_of_bits/(number_of_bits+1); %include leakage power
P_logic_phi = (log2(nb_MAC) + 1) * nb_MAC * 8*C_logic * Vref^2 * clk_freq;  %alpha logic assumed 1 for shift register
P_logic = P_logic_SAR+P_logic_phi;

P_switch = 4*C_logic* (2*nb_activations)* Vref^2*clk_freq;

V_fs = Vref;
C_load_sar = 4*C_logic;  %Assume latch used for input to SAR
E_comp = 2*number_of_bits*log(2)*V_eff*V_fs*C_load_sar;
P_comp = nb_channels/(nb_MAC) * clk_freq * number_of_bits/(number_of_bits+1) * E_comp;

P_tot_SAR = (P_dac+P_comp+P_logic);

P_LNA = Vref*I_LNA;

I_SH = 2*sampling_gain*minimum_sampling_cap*Vref*clk_freq*nb_channels/(nb_MAC);
P_SH = Vref * I_SH;

% Transmission power

transmission_power_per_bit = 10^-9;
P_transmission_full_signal = clk_freq/nb_MAC * nb_channels * number_of_bits * transmission_power_per_bit;


P_tot = P_tot_SAR + P_LNA + P_SH + P_transmission_full_signal + P_switch;
P_LNA_percentage = 100*P_LNA/P_tot;
P_SH_percentage = 100*P_SH/P_tot;
P_comp_percentage = 100*P_comp/P_tot;
P_dac_percentage = 100*P_dac/P_tot;
P_SAR_percentage = 100*P_tot_SAR/P_tot;
P_switch_percentage = 100*P_switch/P_tot;
P_logic_percentage = 100*P_logic/P_tot;
P_transmission_full_percentage = 100*P_transmission_full_signal/P_tot;

%% Calculate total capacitance/area required - trade-off multi-ADC VS OC-ADC
%Total C & overclocking ratio

min_nb_ADC = ceil(nb_activations*(number_of_bits+3)/fmax_ADC);
C_MAC = C_sense_size_multiplier*nb_activations + C_sense_size_multiplier*cap_ratio*nb_channels;
C_SAR = 2^number_of_bits * (min_nb_ADC);
total_C = C_MAC*MAC_cap + C_SAR*DAC_cap + (Cs_unit_cap + Cf_unit_cap)*MAC_cap*min_nb_ADC;
OCR = nb_activations*(number_of_bits+3)/(min_nb_ADC);
        
total_C_unit = round(total_C/minimum_technology_cap);
area_C = total_C*area_per_C;

C_MAC_percentage = 100*(C_MAC*MAC_cap/(minimum_technology_cap*total_C_unit));
C_SH_percentage =  100*((Cs_unit_cap + Cf_unit_cap)*MAC_cap*min_nb_ADC/(minimum_technology_cap*total_C_unit));

C_SAR_percentage = 100*( C_SAR*DAC_cap /(minimum_technology_cap*total_C_unit));

%% Figure plotting

% if (ssim_1_LC>ssim_2)
%     SecondPlot = recovered_signal_LC;
%     mse_SecondPlot = mse_out_LC;
%     ssim_SecondPlot = ssim_1_LC;
% else
%     SecondPlot = recovered_signal_2;
%     mse_SecondPlot = mse_out_2;
%     ssim_SecondPlot = ssim_2;
% end
% Reconstruction
% f = figure;
% f.Position = [300 300 1000 800];
% subplot(221)
% hold on
% plot(recovered_signal)
% plot(total_gain*ADC_input)
% plot(SecondPlot)
% hold off
% legend("Recovered values", "Input values", "Leakage corrected recovery");
% title("Recovered vs input: MSE = " + mse_out + ". SSIM = "+ssim_1)
% subplot(222)
% hold on
% plot(expected_recovered_signal)
% plot(recovered_signal)
% plot(recovered_signal_LC)
% hold off
% legend("Digital recovery", "Analog recovery");
% title("Digital recovery vs analog recovery, MSE = " + mse_out_expected + ". SSIM = "+ssim_expected)
% 
% % Plot expected coeff (previous model) vs actual coeff
% subplot(223)
% hold on
% plot(expected_coeff)
% plot(cs_out)
% plot(cs_out_LC)
% hold off
% title("Digital coefficients vs analog coefficients");
% legend("Digital coefficients", "Analog coefficients");
% subplot(224)
% plot(expected_coeff-cs_out)
% title("Difference between ideal and calculated coefficients")
% 
% % Power & area
% figure
% hold on
% subplot(311)
% bar([P_LNA_percentage; P_SAR_percentage; P_transmission_full_percentage; P_SH_percentage; P_logic_percentage; P_switch_percentage; P_dac_percentage]);
% title("Power per block in percentage of total. Total power = " + P_tot);
% set(gca, 'XTick', 1:7, 'XTickLabels', {"LNA", "SAR", "Transmission", "S&H", "Logic", "Switches", "DAC"});
% subplot(312)
% plot(total_C_unit)
% title("Area per configuration. Chosen config of " +nb_ADC+" ADCs has an approx area of " + area_C(nb_ADC)*10^(12) + "um^2.");
% ylabel("Total amount of unit caps");
% xlabel("Number of ADCs");
% set(gca, 'XTick', min_nb_ADC:min_nb_ADC+3);
% xline(nb_ADC)
% subplot(313)
% plot(OCR)
% set(gca, 'XTick', min_nb_ADC:min_nb_ADC+3);
% xlabel("Number of ADCs");
% ylabel("Overclocking ratio");
% xline(nb_ADC)
% title("Overclocking ratio for given number of ADCs")


%% Plot for thesis
% figure
% hold on
% plot(total_gain*ADC_input)
% plot(recovered_signal_LC)
% legend("Input signal", "Recovered signal")
% title("NMSE = "+num2str(mse_out_LC) + ", SSIM = "+num2str(ssim_1_LC))
% xlabel("Sample index")
% ylabel("Signal (V)")
% hold off
% plot_paper
% 
figure
hold on
plot(expected_coeff)
plot(cs_out)
plot(cs_out_LC)
legend("Ideal coefficients", "Output coefficients", "Leakage corrected coefficients")
xlabel("Coefficient index")
ylabel("Signal (V)")
hold off
plot_paper

figure
hold on
plot(expected_coeff-cs_out_LC)
xlabel("Coefficient index")
ylabel("Coefficient error (V)")
hold off
plot_paper

% figure
% hold on
% plot(expected_recovered_signal)
% plot(recovered_signal_LC)
% legend("Ideal recovery", "Leakage corrected recovery")
% title("Degradation w.r.t NMSE = " + num2str(mse_out_LC-mse_out_expected) + " , SSIM = " + num2str(ssim_expected - ssim_1_LC))
% hold off
% xlabel("Sample index")
% ylabel("Signal (V)")
% plot_paper
% 
% figure
% hold on
% bar([P_SH_percentage; P_dac_percentage;P_logic_percentage; P_comp_percentage;P_switch_percentage ]);
% set(gca, 'XTick', 1:5, 'XTickLabels', { "S&H", "DAC", "Logic", "Comparator", "Switching network"});
% ylabel("Relative power draw percentages without Tx & LNA")
% title("Total power: " + P_tot*1e6 + "\muW")
% plot_paper
% axes('Position',[.7 .7 .2 .2])
% box on
% bar([P_LNA_percentage; 100-(P_LNA_percentage+P_transmission_full_percentage); P_transmission_full_percentage; ]);
% set(gca, 'XTick', 1:3, 'XTickLabels', {"LNA", "Other", "Transmission"});
% title("With TX and LNA")
% plot_paper

% figure
% hold on
% bar(1/100*[C_SH_percentage*total_C_unit; C_MAC_percentage*total_C_unit;C_SAR_percentage*total_C_unit]);
% set(gca, 'XTick', 1:3, 'XTickLabels', { "S&H", "MAC array", "SAR"});
% ylabel("Number of unit capacitors")
% % title("Total unit capacitance: " + total_C_unit)
% plot_paper
% 


% figure
% hold on
% plot(expected_coeff, "--")
% scatter(1:nb_channels, cs_out, 50, "filled")
% scatter(1:nb_channels, cs_out_LC, 50, "LineWidth", 2)
% xlabel("Coefficient index")
% ylabel("Signal (V)")
% hold off
% plot_paper