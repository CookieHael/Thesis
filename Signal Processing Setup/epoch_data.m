function [nb_epochs,epoched_data] = epoch_data(data, nb_MAC)


if rem(length(data), nb_MAC) > 0
    data = [data, zeros(1, nb_MAC-rem(length(data), nb_MAC))];
end
nb_epochs = length(data)/nb_MAC;
epoched_data = reshape(data, nb_MAC, nb_epochs);

end

