function [rms, nmse] = calculateRMS(input_values, output_values)
    if (length(input_values) ~= length(output_values))
        rms = -inf;
        mse = -inf;
        return
    end

    rms = sqrt(  mean( (input_values-output_values).^2 ) );
    nmse = (norm(input_values - output_values,'fro')/norm(input_values,'fro'))^2;
end