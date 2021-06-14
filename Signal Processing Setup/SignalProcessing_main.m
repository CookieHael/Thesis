% clear;
% verbose = false;
%  tic
% 
% %%  Constants 
% k = 1.38064852e-23;
% T = 293;
% V_thermal = 25.27*10^-3;
% minimum_technology_cap = 1.995 * 1e-15;
% V_eff = 0.1;
% area_per_C = 0.001025;  %F/m^2
% 
% %% Process parameters
% Vref = 1;
% gmoverid = 20; %%I=150nA
% fmax_ADC = 10e6;   %% Determine from simulation!
% 
% %% System settings
% number_of_bits = 8;
% input_BW = 128; %This should be the Nyquist frequency, and input should be oversampled wrt this
% signal_peak_amplitude = .001;
% signal_bias = 0;
% nb_MAC = 512;
% cap_ratio = 15;
% nb_channels = 256;
% nb_activations = 2;
% 
% non_idealities_on = 1;
% noise_on = 1;
% 
% if ~noise_on
%     warning("Noise currently off!")
% end
% 
% noise_quantization = (Vref)^2/12*2^(-2*number_of_bits);
% 
% %% Sensing matrix setup
% 
% if(nb_channels == 150 && nb_MAC == 384)
%     sensing_matrix_large = load("sensing150.mat").sensing_matrix_large;  %%Capacitor ratio was 35
%     sensing_matrix_corrected = load("sensing150.mat").sensing_matrix_corrected;
% elseif (nb_channels == 192 && nb_MAC == 384)
%     sensing_matrix_large = load("sensing192.mat").Phi;
%     sensing_matrix_corrected = load("sensing192.mat").corrected_Phi;
% else
%     [sensing_matrix_large, sensing_matrix_corrected] = generateCorrectedSRBM(nb_channels, nb_MAC, nb_activations, cap_ratio);
% end
% if (cap_ratio ~=35)
%     sensing_matrix_corrected = rescaleSRBM(sensing_matrix_large, cap_ratio);
% end
% 
% if( nb_MAC ~= size(sensing_matrix_large, 2) || nb_channels ~= size(sensing_matrix_large,1))
%     error("Wrong matrix dimensions");
% end
% 
% sensing_order = determineSensingOrder(sensing_matrix_large, nb_activations, 1);
% sensing_matrix_sensing_order = sensing_matrix_large(sensing_order,:);
% sensing_matrix_prep = sensing_matrix_sensing_order';
% sensing_matrix = sensing_matrix_prep(:)';
% 
% %% Clock & sine timing settings
% 
% clk_freq = 2.1*input_BW;
% clk_period = 1/clk_freq;
% 
% %% Sample time jitter
% sampling_jitter_on = noise_on;
% sample_time_jitter = 16e-12;
% 
% %% Capacitor sizing
% minimum_sampling_cap = 12*k*T*2/3*2^(2*number_of_bits)/Vref^2;
% minimum_cap = max(minimum_technology_cap,minimum_sampling_cap);
% 
% %% MAC non-idealities & noise
% C_sense_size = 1;
% C1 = C_sense_size*minimum_cap;
% C2 = C1*cap_ratio; 
% MM_incomplete_transfer_coeff = C1/(C1+C2);
% if (non_idealities_on)
%     MM_attenuation_coeff = C2/(C1+C2);  %Badly chosen name: this is voltage attenuation by charge sharing architecture
% else
%     MM_attenuation_coeff = 1;  %Actually more like 1-leakage
% end
% 
% MM_sampling_noise = noise_on*((k*T * C1/(C1+C2)^2) + (k*T*(C1*C2/(C1+C2))/C2^2));  % Formula with square root is sigma, model uses variance
% 
% %% Sample & hold (with gain)
% switch_noise_on = noise_on;
% Cs_unit_cap = 20/MM_incomplete_transfer_coeff;
% Cf_unit_cap = 20;
% sampling_gain = Cs_unit_cap/Cf_unit_cap;
% switching_noise = k*T/(Cf_unit_cap*minimum_cap);
% 
% %% Comparator offset/hysteresis & noise
% comparator_offset_on = noise_on;
% comparator_high_offset = 0.2*10e-4;
% comparator_low_offset = -.1*10e-4;
% comparator_noise_on = noise_on;
% comparator_noise_rms = 1*k*T* 2/3 / minimum_cap;
% 
% %% LNA
% LNA_noise_on = noise_on;
% LNA_gain = 100;
% LNA_bandwidth = 3*input_BW; %% in Hz
% LNA_NEF = 1.08;
% 
% % Determine minimum current if bandwidth-constrained
% gm_min_1 = LNA_bandwidth*2*pi*nb_activations*minimum_cap;
% Id_min_1 = gm_min_1*gmoverid;
% LNA_input_referred_noise_1 = LNA_NEF/(sqrt(2*Id_min_1/(pi*V_thermal*4*k*T*LNA_bandwidth)));
% 
% % Determine minimum current if slewrate-constrained
% SR_required = Vref*clk_freq; %Minimum current
% Id_min_2 = SR_required*nb_activations*minimum_cap;
% gm_min_2 = Id_min_2*gmoverid;
% LNA_input_referred_noise_2 = LNA_NEF/(sqrt(2*Id_min_2/(pi*V_thermal*4*k*T*LNA_bandwidth)));
% 
% % Determine minimum current if noise-constrained
% LNA_input_referred_noise_3 = noise_quantization*10;
% Id_min_3 = (LNA_NEF/LNA_input_referred_noise_3)^2*pi*4*k*T*LNA_bandwidth*V_thermal;
% 
% % Determine which of the three cases limits LNA performance
% I_LNA = max([Id_min_1, Id_min_2, Id_min_3]);
% if (I_LNA==Id_min_3 && Id_min_3>130e-09)
%     warning("gm/Id not really valid anymore, 20 at 130nA but lower at higher currents");
% end
% 
% LNA_in_rms = min([LNA_input_referred_noise_1,LNA_input_referred_noise_2,LNA_input_referred_noise_3]); % Take maximum current of 3 required, noise will be minimum due to NEF formula
% LNA_SR = I_LNA/(nb_activations*minimum_cap);
% pd = makedist('Normal', 0, LNA_in_rms);
% rng(42069);
% LNA_noise_vector = [pd.random(nb_MAC,1);0]';
% 
% %% Total gain check
% total_gain = LNA_gain*sampling_gain*MM_incomplete_transfer_coeff;
% if(signal_peak_amplitude*total_gain>Vref/2)
%     warning("Careful! Input signal might be saturating ADC conversion range.\n");
% end
% 
% %% ZOH settings
% 
% ZOH_delta = 0*Vref/(2^number_of_bits);
% if (ZOH_delta>0)
%     fprintf("Reminder: ZOH is on, signal quality may be degraded");
% end
% %% Inputs
% names = ['A', 'B', 'C', 'D', 'E'];
% names2 = ['Z', 'O', 'N', 'F', 'S'];
% 
% %Set up model for parsim
% model = "MMSAR_Sequential";
% load_system(model);
% for index=1:5
%     index
%     tic
%     simIn = Simulink.SimulationInput(model);
% 
%     data = load_EEG_data(names(index));
%     ADC_input_reshaped = form_input(data, nb_MAC, nb_channels);
%     sim_number = index;
%     [nb_epochs, epoched_data] = epoch_data(data, nb_MAC);
% 
%     sim_out=sim("MMSAR_Sequential",(nb_channels+1)*clk_period*nb_MAC*nb_epochs);
% 
%     [range, minima] = getRange(names(index));
%     %% Processing
%     output = zeros(nb_epochs, nb_MAC);
%     rms_out = zeros(5, nb_epochs);
%     nmse_out = zeros(5, nb_epochs);
% 
%     % Simulate
%     
%     for j = 1:nb_epochs
%         cs_out = zeros(nb_channels,1);
%         for i = 1:nb_channels
%             cs_out(sensing_order(i)) = sim_out(1,index).yout.signals(1).values(i+3+(j-1)*nb_channels);
%         end
% 
%         %nb_ZOH = sim_out.yout.signals(2).values(end);
% 
%         % Compressive sensing reconstruction      
%         if (non_idealities_on)
%             A = sensing_matrix_corrected*wmpdictionary(nb_MAC, 'lstcpt', {'dct'});
%         else
%             A = sensing_matrix_large*wmpdictionary(nb_MAC, 'lstcpt', {'dct'});
%         end
% 
%         recovery = BSBL_BO(A, cs_out, 1:15:nb_MAC, 0, 'prune_gamma',-1, 'max_iters',10);
%         recovered_signal = (idct(recovery.x))';
%         recovered_signal = recovered_signal(1:nb_MAC);
% 
%         [rms_out_partial, nmse_out_partial] = calculateRMS(epoched_data(:,j)'*total_gain, recovered_signal);
%        
%         %% Output saving
%         
%         rms_out(index, j) = rms_out_partial;
%         nmse_out(index, j) = nmse_out_partial;
%         %Restoring to original range
% 
%         recovered_signal = 10*recovered_signal;
%         recovered_signal = (recovered_signal+1)/2;
%         recovered_signal = recovered_signal*range(ceil(j/8)) + minima(ceil(j/8));
%         output(j,:) = recovered_signal;
%         
%     end
%     
%     for file_index = 1:100
% 
%         if (file_index<10)
%             fid = fopen("Output data/" + names(index) + "/" + names2(index) + "00" + int2str(file_index)+".txt", 'w');
%         else
%             fid = fopen("Output data/" + names(index) + "/" + names2(index) + "0" + int2str(file_index)+".txt", 'w');
%         end
% 
%         for j = 1:8
%         fprintf(fid, '%f \n', output(8*(file_index-1)+j,:));
%         end
% 
%     end
% toc
% end
% fclose("all");
% toc



clear;
verbose = false;
tic

%%  Constants 
k = 1.38064852e-23;
T = 293;
V_thermal = 25.27*10^-3;
minimum_technology_cap = 1.995 * 1e-15;
V_eff = 0.1;
area_per_C = 0.001025;  %F/m^2

%% Process parameters
Vref = 2;
gmoverid = 20; %%I=150nA
fmax_ADC = 10e6;   %% Determine from simulation!

%% System settings
number_of_bits = 8;
input_BW = 128; %This should be the Nyquist frequency, and input should be oversampled wrt this
signal_peak_amplitude = .001;
signal_bias = 0;
nb_MAC = 512;
cap_ratio = 15;
nb_channels = 256;
nb_activations = 2;

non_idealities_on = 1;
noise_on = 1;

if ~noise_on
    warning("Noise currently off!")
end

noise_quantization = (Vref)^2/12*2^(-2*number_of_bits);

%% Sensing matrix setup

if(nb_channels == 150 && nb_MAC == 384)
    sensing_matrix_large = load("sensing150.mat").sensing_matrix_large;  %%Capacitor ratio was 35
    sensing_matrix_corrected = load("sensing150.mat").sensing_matrix_corrected;
elseif (nb_channels == 192 && nb_MAC == 384)
    sensing_matrix_large = load("sensing192.mat").Phi;
    sensing_matrix_corrected = load("sensing192.mat").corrected_Phi;
else
    [sensing_matrix_large, sensing_matrix_corrected] = generateCorrectedSRBM(nb_channels, nb_MAC, nb_activations, cap_ratio);
end
if (cap_ratio ~=35)
    sensing_matrix_corrected = rescaleSRBM(sensing_matrix_large, cap_ratio);
end

if( nb_MAC ~= size(sensing_matrix_large, 2) || nb_channels ~= size(sensing_matrix_large,1))
    error("Wrong matrix dimensions");
end

sensing_order = determineSensingOrder(sensing_matrix_large, nb_activations, 1);
sensing_matrix_sensing_order = sensing_matrix_large(sensing_order,:);
sensing_matrix_prep = sensing_matrix_sensing_order';
sensing_matrix = sensing_matrix_prep(:)';
clear sensing_matrix_prep

%% Clock & sine timing settings

clk_freq = 2.1*input_BW;
clk_period = 1/clk_freq;

%% Sample time jitter
sampling_jitter_on = noise_on;
sample_time_jitter = 16e-12;

%% Capacitor sizing
minimum_sampling_cap = 12*k*T*2/3*2^(2*number_of_bits)/Vref^2;
minimum_cap = max(minimum_technology_cap,minimum_sampling_cap);

%% MAC non-idealities & noise
C_sense_size = 1;
C1 = C_sense_size*minimum_cap;
C2 = C1*cap_ratio; 
MM_incomplete_transfer_coeff = C1/(C1+C2);
if (non_idealities_on)
    MM_attenuation_coeff = C2/(C1+C2);  %Badly chosen name: this is voltage attenuation by charge sharing architecture
else
    MM_attenuation_coeff = 1;  %Actually more like 1-leakage
end

MM_sampling_noise = noise_on*((k*T * C1/(C1+C2)^2) + (k*T*(C1*C2/(C1+C2))/C2^2));  % Formula with square root is sigma, model uses variance

%% Sample & hold (with gain)
switch_noise_on = noise_on;
Cs_unit_cap = 20/MM_incomplete_transfer_coeff;
Cf_unit_cap = 20;
sampling_gain = Cs_unit_cap/Cf_unit_cap;
switching_noise = k*T/(Cf_unit_cap*minimum_cap);

%% Comparator offset/hysteresis & noise
comparator_offset_on = noise_on;
comparator_high_offset = 0.2*10e-4;
comparator_low_offset = -.1*10e-4;
comparator_noise_on = noise_on;
comparator_noise_rms = 1*k*T* 2/3 / minimum_cap;

%% LNA
LNA_noise_on = noise_on;
LNA_gain = 100;
LNA_bandwidth = 3*input_BW; %% in Hz
LNA_NEF = 1.08;

% Determine minimum current if bandwidth-constrained
gm_min_1 = LNA_bandwidth*2*pi*nb_activations*minimum_cap;
Id_min_1 = gm_min_1*gmoverid;
LNA_input_referred_noise_1 = LNA_NEF/(sqrt(2*Id_min_1/(pi*V_thermal*4*k*T*LNA_bandwidth)));

% Determine minimum current if slewrate-constrained
SR_required = Vref*clk_freq; %Minimum current
Id_min_2 = SR_required*nb_activations*minimum_cap;
gm_min_2 = Id_min_2*gmoverid;
LNA_input_referred_noise_2 = LNA_NEF/(sqrt(2*Id_min_2/(pi*V_thermal*4*k*T*LNA_bandwidth)));

% Determine minimum current if noise-constrained
LNA_input_referred_noise_3 = noise_quantization*10;
Id_min_3 = (LNA_NEF/LNA_input_referred_noise_3)^2*pi*4*k*T*LNA_bandwidth*V_thermal;

% Determine which of the three cases limits LNA performance
I_LNA = max([Id_min_1, Id_min_2, Id_min_3]);
if (I_LNA==Id_min_3 && Id_min_3>130e-09)
    warning("gm/Id not really valid anymore, 20 at 130nA but lower at higher currents");
end

LNA_in_rms = min([LNA_input_referred_noise_1,LNA_input_referred_noise_2,LNA_input_referred_noise_3]); % Take maximum current of 3 required, noise will be minimum due to NEF formula
LNA_SR = I_LNA/(nb_activations*minimum_cap);
pd = makedist('Normal', 0, LNA_in_rms);
rng(42069);
LNA_noise_vector = [pd.random(nb_MAC,1);0]';

%% Total gain check
total_gain = LNA_gain*sampling_gain*MM_incomplete_transfer_coeff;
if(signal_peak_amplitude*total_gain>Vref/2)
    warning("Careful! Input signal might be saturating ADC conversion range.\n");
end

%% ZOH settings

ZOH_delta = 0*Vref/(2^number_of_bits);
if (ZOH_delta>0)
    fprintf("Reminder: ZOH is on, signal quality may be degraded");
end
%% Inputs
names = ['A', 'B', 'C', 'D', 'E'];
names2 = ['Z', 'O', 'N', 'F', 'S'];

%Set up model for parsim
model = "MMSAR_Sequential2019";
load_system(model);
save('simVariables');

simIn(1:500) = Simulink.SimulationInput(model);
for i=1:5
    for j=1:100
        data = load_EEG_data(names(i), j);
        [nb_epochs, epoched_data] = epoch_data(data, nb_MAC);
        simIn((i-1)*100+j) = simIn((i-1)*100+j).setVariable("epoched_data", epoched_data);
        simIn((i-1)*100+j) = simIn((i-1)*100+j).setVariable("sim_number", (i-1)*100+j);
        simIn((i-1)*100+j) = simIn((i-1)*100+j).setVariable("nb_epochs", nb_epochs);
        simIn((i-1)*100+j) = simIn((i-1)*100+j).setModelParameter("StartTime", "0", "StopTime", num2str((nb_channels+1)*clk_period*nb_MAC*nb_epochs));
    end
end
toc
sim_out=parsim(simIn, 'showProgress', 'on', 'SetupFcn', @() (evalin('base', 'load simVariables.mat')));
save('simOut_'+num2str(nb_channels)+'ch.mat', 'sim_out');
toc

parfor index=1:5
    index
    for l = 1:100
    l

        data = load_EEG_data(names(index), l);
        [range, minima] = getRange(names(index), l);
        [nb_epochs, epoched_data] = epoch_data(data, nb_MAC);
        %% Processing
        output = zeros(800, nb_MAC);
        rms_out = zeros(5, 800);
        nmse_out = zeros(5, 800);

        if (l<10)
            fid = fopen("Output data/" + names(index) + "/" + names2(index) + "00" + int2str(l)+".txt", 'w');
        elseif (l==100)
            fid = fopen("Output data/" + names(index) + "/" + names2(index) + "0" + int2str(l)+".txt", 'w');
        else
            fid = fopen("Output data/" + names(index) + "/" + names2(index) + "0" + int2str(l)+".txt", 'w');
        end
        
        % Simulate
        for epoch=1:nb_epochs
            cs_out = zeros(nb_channels,1);
            for i = 1:nb_channels
                cs_out(sensing_order(i)) = sim_out(1,100*(index-1)+l).yout.signals(1).values(i+3+nb_channels*(epoch-1));
            end
           
            %nb_ZOH = sim_out.yout.signals(2).values(end);

            % Compressive sensing reconstruction      
            if (non_idealities_on)
                A = sensing_matrix_corrected*wmpdictionary(nb_MAC, 'lstcpt', {'dct'});
            else
                A = sensing_matrix_large*wmpdictionary(nb_MAC, 'lstcpt', {'dct'});
            end

            recovery = BSBL_BO(A, cs_out, 1:15:nb_MAC, 0, 'prune_gamma',-1, 'max_iters',10);
            recovered_signal = (idct(recovery.x))';
            recovered_signal = recovered_signal(1:nb_MAC);

            [rms_out_partial, nmse_out_partial] = calculateRMS(epoched_data(:,epoch)'*total_gain, recovered_signal);
            nmse_out_partial

            %% Output saving

            rms_out(index, epoch+(l-1)*8) = rms_out_partial;
            nmse_out(index, epoch+(l-1)*8) = nmse_out_partial;
            %Restoring to original range and writing

            recovered_signal = 10*recovered_signal;
            recovered_signal = (recovered_signal)/2;
            recovered_signal = recovered_signal*range + minima;
            output(epoch+(l-1)*8,:) = round(recovered_signal);
            
            fprintf(fid, '%f \n', recovered_signal);
            
        end 

    end
toc
end
fclose("all");
toc

