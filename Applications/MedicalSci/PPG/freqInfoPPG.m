function [posl,posh,shat] = freqInfoPPG(sig)
% posl for lung
% posh for heart

N = length(sig);
if mod(N,2) > 0.5
    sig = sig(1:end-1);
    N = N -1;
end
shat = fft(sig);
shat = abs(shat(1:N/2));
shat(1:10) = 0;
[~,pos1] = max(shat);
shat(max(1,round(pos1*0.8)):round(pos1*1.2)) = 0;
[~,pos2] = max(shat);
if pos1(1) < pos2(1)
    temp = pos1(1);
    pos1 = pos2(1);
    pos2 = temp;
end
if abs(pos1(1)/pos2(1)-2) < 0.5
    posh = pos2(1);
    shat(round(pos2(1)*0.8):end) = 0;
    [~,posl] = max(shat);
else
    posh = pos1(1);
    posl = pos2(1);
end
end



