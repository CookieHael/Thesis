function predictedX = findMismatch(input, sensing_vector, output_coeff, gain)

idx = find(sensing_vector);
p = zeros(length(idx)+1, 1);
p(1) = 1;
p(end) = (output_coeff/(input*gain)-1);
predictedX = roots(p');


end

