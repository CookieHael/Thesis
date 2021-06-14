load 'simVariables';
load 'simOut';
tic;
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

            %% Output saving

            rms_out(index, epoch+(l-1)*8) = rms_out_partial;
            nmse_out(index, epoch+(l-1)*8) = nmse_out_partial;
            %Restoring to original range and writing

            recovered_signal = 10*recovered_signal;
            recovered_signal = (recovered_signal)/2;
            recovered_signal = recovered_signal*range + minima;
            
            fprintf(fid, '%f \n', recovered_signal);
            
        end 

    end

end
fclose("all");
toc;