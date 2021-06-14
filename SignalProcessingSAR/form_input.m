function output = form_input(data,nb_MAC, nb_channels)

if rem(length(data), nb_MAC) > 0
    data = [data, zeros(1, nb_MAC-rem(length(data), nb_MAC))];
end
output = reshape(data, nb_MAC, length(data)/nb_MAC);
output = repelem(output, 1, nb_channels);
output = output(:);
end

