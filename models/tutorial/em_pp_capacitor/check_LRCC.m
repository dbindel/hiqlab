% -- Assume the following exists in memory
%
%   -Admittance of model [freq,H] 
%   -Equivalent LRCC parmaters [eL,eR,eC,eC0]

Hl   = zeros(1,length(freq));
for j=1:length(freq)
    w    = freq(j);
    Hl(j)= i*w*eC0 + 1/(i*w*eL + eR + 1/(i*w*eC));
end;

figure(1)
plot(freq/2/pi,20*log10(abs(H)),'r-');
hold on;
plot(freq/2/pi,20*log10(abs(Hl)),'bd');
hold off;

figure(2);
plot(freq/2/pi,angle(H),'r-');
hold on;
plot(freq/2/pi,angle(Hl),'bd');
hold off;
