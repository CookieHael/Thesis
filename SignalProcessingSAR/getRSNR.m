function [RSNR] = getRSNR(original,recovery)
%GETRSNR Calculate reconstruction SNR

   RSNR=20*log10(norm(original)/norm(original-recovery));


end

