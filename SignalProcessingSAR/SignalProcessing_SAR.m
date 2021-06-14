clear;

verbose = false;

%%  ----------  Circuit parameters  --------------
tic
for nb_bits=4:12
    number_of_bits = nb_bits
    input_BW=256;
    V_thermal = 25.27*10^-3;
    k = 1.38064852e-23;
    F_per_M = 0.001025;
    T = 293;
    sine_amplitude = .001;
    sine_bias = 0;
    Vref = 2;
    noise_quantization = sqrt(Vref^2/12*2^(-2*number_of_bits));
    C_mismatch_parameter = 3.4878e-09;


    % Clock & sine timing settings
    sine_frequency_rad = 2*pi*input_BW;  % Set to 256Hz
    %nyquist_frequency = sine_frequency*number_of_bits*2;
    clk_freq = 2.1*(number_of_bits+1)*sine_frequency_rad;
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
    comparator_high_offset = 8e-4;
    comparator_low_offset = 8e-4;
    comparator_noise_on = non_idealities_on;
    comparator_noise_rms = 1*k*T* 2/3 / minimum_cap;



    %% LNA
    LNA_noise_on = non_idealities_on;
    LNA_gain = 100;
    LNA_bandwidth = 30*input_BW; %% in Hz
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
    total_gain = sampling_gain*LNA_gain;
    
    if(total_gain*sine_amplitude>Vref)
        warning("Careful! Input signal might be saturating ADC conversion range.\n");
    end

    %% Inputs
    names = ['A', 'B', 'C', 'D', 'E'];
    names2 = ['Z', 'O', 'N', 'F', 'S'];

    %Set up model for parsim
    model = "ADC2019";
    load_system(model);
    save('simVariables');

    simIn(1:500) = Simulink.SimulationInput(model);
    for i=1:5
        for j=1:100
            data = load_EEG_data(names(i), j);
            simIn((i-1)*100+j) = simIn((i-1)*100+j).setVariable("ADC_input", data);
            simIn((i-1)*100+j) = simIn((i-1)*100+j).setVariable("sim_number", (i-1)*100+j);
            simIn((i-1)*100+j) = simIn((i-1)*100+j).setModelParameter("StartTime", "0", "StopTime", num2str((length(data)+1)*(number_of_bits+1)*clk_period));
        end
    end
    toc
    sim_out=parsim(simIn, 'showProgress', 'on', 'SetupFcn', @() (evalin('base', 'load simVariables.mat')));
    save("simOutSAR_"+int2str(number_of_bits)+"ch.mat", 'sim_out');
    toc

    for index=1:5
        index
        for l = 1:100
        

            data = load_EEG_data(names(index), l);
            [range, minima] = getRange(names(index), l);
            %% Processing
            output = zeros(500, length(data));
            rms_out = zeros(5, 500);
            nmse_out = zeros(5, 500);
            RSNR_out = zeros(5, 500);
            ssim_out = zeros(5,500);
            
            if (~exist("Output data"+int2str(number_of_bits)+"/" + names(index), 'dir'))
                mkdir("Output data"+int2str(number_of_bits)+"/" + names(index), 'dir')
            end
            
            if (l<10)
                fid = fopen("Output data"+int2str(number_of_bits)+"/" + names(index) + "00" + int2str(l) + "SAR.txt", 'w');
            elseif (l==100)
                fid = fopen("Output data"+int2str(number_of_bits)+"/" + names(index) + "0" + int2str(l) + "SAR.txt", 'w');
            else
                fid = fopen("Output data"+int2str(number_of_bits)+"/" + names(index) + "0" + int2str(l) + "SAR.txt", 'w');
            end

            recovered_signal = sim_out(1,100*(index-1)+l).yout.signals(1).values(3:end);


            [rms, nmse] = calculateRMS(data'*total_gain, recovered_signal);
            RSNR = getRSNR(data'*total_gain, recovered_signal);
            windowLen = 100;
            [mssim, ssim_map] = ssim_1d(data'*total_gain, recovered_signal, windowLen);


            %% Output saving

            rms_out(index, l) = rms;
            nmse_out(index, l) = nmse;
            RSNR_out(index, l) = RSNR;
            ssim_out(index, l) = mssim;
            %Restoring to original range and writing

            recovered_signal = 10*recovered_signal;
            recovered_signal = (recovered_signal)/2;
            recovered_signal = recovered_signal*range + minima;
            output(l+(index-1)*100,:) = round(recovered_signal);

            fprintf(fid, '%f \n', recovered_signal);


        end
    toc
    end
    fclose("all");
    toc
end
