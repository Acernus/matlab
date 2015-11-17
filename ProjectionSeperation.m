clc;
clear all;
B1 = [0 : 0.001 : 10];
B2 = [0 : 0.001 : 10];

low_spec = load('80kv_Spectrum.dat');
high_spec = load('150kv_Spectrum.dat');
%80kv
lu1 = 2.037 * 10^-2;% Carbon    cm^2/g
lu2 = 1.104 * 10^-1;% Iron      cm^2/g
%150kv
hu1 = 2.449 * 10^-2;
hu2 = 7.961 * 10^-2;
[h1, w1] = size(B1);
[h2, w2] = size(B2);
[hs1, ws1] = size(low_spec);
[hs2, ws2] = size(high_spec);
PL = [];
PH = [];
for i = 1 : w1
    for j = 1 : w2
        %calculate value of PL. 
        lowlength = 0;
        lowlengthMultiple = 0;
        for m = 1 : hs1
            lowlength = lowlength + low_spec(m, 2);
            lowlengthMultiple = lowlengthMultiple + low_spec(m, 2) * exp(-B1(i) * lu1 - B2(j) * lu2);
        end
        PL = [PL, -log(lowlengthMultiple) + log(lowlength)];
        %calculate value of PH
        highlength = 0;
        highlengthMultiple = 0;
        for m = 1 : hs2
            highlength = highlength + high_spec(m, 2);
            highlengthMultiple = highlengthMultiple + high_spec(m, 2) * exp(-B1(i) * hu1 - B2(j) * hu2);
        end
        PH = [PH, -log(highlengthMultiple) + log(highlength)];
    end
end
