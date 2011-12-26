% -- Assume the following exists in memory
%
%   -Transfer of model [freq,H], Rs,Rf

% -- Equivalent LRCC parameters of resonator 1
eL1 = 1.793695949424715e-01; 
eR1 = 5.368942077626494e+02;
eC1 = 7.202930365547299e-18;
eC01= 3.408078874904376e-13;

% -- Equivalent LRCC parameters of resonator 2
eL2 = 1.793732685597236e-01;
eR2 = 5.368998301008118e+02;
eC2 = 7.202926816733675e-18;
eC02= 3.408078874904376e-13;

% -- Resistances
Rf  = 1000;
Rs  = 10000;

Hl   = zeros(1,length(freq));
H1   = zeros(1,length(freq));
H2   = zeros(1,length(freq));

for j=1:length(freq)

    w    = freq(j)*2*pi;

    Y1   = i*w*eC01 + 1/(i*w*eL1 + eR1 + 1/(i*w*eC1));
    Y2   = i*w*eC02 + 1/(i*w*eL2 + eR2 + 1/(i*w*eC2));
    Y3   = 1/Rs + Y2;

    Zt   = Rf + 1/Y1 + 1/Y3;
    Zt1  = Rf + 1/Y1 + Rs;
    Zt2  = Rf        + 1/Y3;

    Hl(j)= (Rs+Rf)/Rs*(1/Y3)/Zt;
    H1(j)= (Rs+Rf)/Rs*    Rs/Zt1;
    H2(j)= (Rs+Rf)/Rs*(1/Y3)/Zt2;
end;

figure(1)
plot(freq/1e6,20*log10(abs(H)),'k-');
hold on;
plot(freq/1e6,20*log10(abs(Hl)),'ro');
plot(freq/1e6,20*log10(abs(H1)),'bd');
plot(freq/1e6,20*log10(abs(H2)),'gd');
hold off;

figure(2);
plot(freq/1e6,angle(H) ,'k-');
hold on;
plot(freq/1e6,angle(Hl),'ro');
plot(freq/1e6,angle(H1),'bd');
plot(freq/1e6,angle(H2),'gd');
hold off;
