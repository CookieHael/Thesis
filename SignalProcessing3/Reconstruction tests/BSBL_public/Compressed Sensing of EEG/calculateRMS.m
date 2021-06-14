function rms_out = calculateRMS(input_values, output_values)
    if (length(input_values) ~= length(output_values))
        rms_out = -inf;
        return
    end
    input_values = input_values(2:length(input_values));
    output_values = output_values(1:length(output_values)-1);
    rms_out = sqrt(  sum( (input_values-output_values).^2 )  / length(output_values) );
    %rms_out = (norm(input_values - output_values,'fro')/norm(input_values,'fro'))^2;
end