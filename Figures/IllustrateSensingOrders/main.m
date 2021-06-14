clear;
leakNew = zeros(100,5);
N = [50, 75, 100, 150, 250];
parfor k = 1:5
    for j=1:1000

        phi = generateSRBM(N(k),384, 2);

        leak_cycles_array=zeros(size(phi,1),1);
        [sensing_matrix, order] = determineSensingOrderNew(phi, 2, 1);

        for i=1:size(phi,1)
                idx = find(phi(i, :));
                idx2 = find(sensing_matrix(i,:));
                if(isempty(idx))
                    idx=384;
                end
                if (isempty(idx2))
                    idx2=384;
                end

                leak_cycles_array(i) =  idx2(1)-idx(1);
        end
        leakNew(j,k) = mean(leak_cycles_array);

    end
end
figure
plot(leakNew)
l=legend("50", "75", "100", "150", "250");
title(l, "M")
plot_paper
xlabel("Sample")
ylabel("Average leakage cycles")

figure
leakNew2=zeros(31,5);
temp = zeros(500,1);
m=1;
for act = 1:5
for k = 50:5:200
    parfor j=1:500

        phi = generateSRBM(k,384, act);

        leak_cycles_array=zeros(size(phi,1),1);
        [sensing_matrix, order] = determineSensingOrderNew(phi, act, 1);

        for i=1:size(phi,1)
                idx = find(phi(i, :));
                idx2 = find(sensing_matrix(i,:));
                if(isempty(idx))
                    idx=384;
                end
                if (isempty(idx2))
                    idx2=384;
                end

                leak_cycles_array(i) =  idx2(1)-idx(1);
        end
        temp(j) = mean(leak_cycles_array);

    end
    leakNew2(m, act) = mean(temp);
    m=m+1;
end
m=1;
end

plot(50:5:200, leakNew2)
plot_paper
l=legend("1", "2", "3", "4", "5");
xlabel("M")
ylabel("Average leakage cycles")
title(l, "#Activations")