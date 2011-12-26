% -- Assume the following exists in memory
%
%   -Transfer of model [freq,H], Rs,Rf
%   -Equivalent LRCC parmaters [eL,eR,eC,eC0]

eL = 1.793695949424714e-01;
eR = 5.368942069576071e+02;
eC = 7.202930365547330e-18;
eC0= 3.408078874904376e-13;
Rs = 50;
Rf = 50;

Hl   = zeros(1,length(freq));
for j=1:length(freq)
    w    = freq(j)*2*pi;
    Y    = i*w*eC0 + 1/(i*w*eL + eR + 1/(i*w*eC));
    Hl(j)= (Rs+Rf)/(Rs+Rf+1/Y);
end;

figure(1)
plot(freq,20*log10(abs(H)),'r-');
hold on;
plot(freq,20*log10(abs(Hl)),'bd');
hold off;

figure(2);
plot(freq,angle(H),'r-');
hold on;
plot(freq,angle(Hl),'bd');
hold off;
