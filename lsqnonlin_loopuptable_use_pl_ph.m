clc;
clear all;
%%
%�������ļ��ж����ߵ����״���SL��SH��
SL = load('newSpectrum80kV.dat');%1-80KeV
SH = load('newSpectrum160kV.dat');%1-160Kev
points = 10;
dLE = 80/points;
dHE = 160/points;
[lht, lwd] = size(SL);
[hht, hwd] = size(SH);
%%
%�ֶβ������˴���ȡ600����
lsamplelen = points;
hsamplelen = points;
sampleSL = [];
sampleSH = [];
for i = 0 : points -1 
    sampleSL = [sampleSL SL(1 + i*dLE)];
    sampleSH = [sampleSH SH(1 + i*dHE)];
end
%������ҵ��ĸߵ������£��������ʵ�˥��ϵ��, size == ����������
%����������ϣ��ֱ�õ��ߵ��������ʵ�˥������
L_Al = interpolant('0_80kv_Al.txt', 3);
L_Carbon = interpolant('0_80kv_Carbon.txt', 3);
H_Al = interpolant('0_160kv_Al.txt', 3);
H_Carbon = interpolant('0_160kv_Carbon.txt', 3);
ul = [];
uh = [];
for i = 0 : lwd
    ulAl = getAttenuationValue(L_Al, 1 + dLE*i);
    ulCarbon = getAttenuationValue(L_Carbon, 1 + dLE*i);
    ul = [ul;ulAl ulCarbon];
end
for i = 0 : hwd
    uhAl = getAttenuationValue(H_Al, 1 + dLE*i);
    uhCarbon = getAttenuationValue(H_Carbon, 1 + dLE*i);
    uh = [uh; uhAl uhCarbon];
end

SLintergration = 0;
SHintergration = 0;

for i = 1 : lsamplelen
    SLintergration = SLintergration + sampleSL(i) * dLE;
end
for i = 1 : hsamplelen
    SHintergration = SHintergration + sampleSH(i) * dHE;
end
%%

%ȡ�������ߵ��ܵ�ͶӰ����
highprj = readbin('PROJ_MATCH_H2.BIN');
lowprj = readbin('PROJ_MATCH_L2.BIN');
%����ͶӰ����ƫ��
for i = 1 : 3
    highprj = [highprj(2:400,:);highprj(1,:)];
end

[height, width] = size(lowprj);
%%
%�ȶ�����������ͶӰͼ�ı���
ul1 = ul(:, 1);
ul2 = ul(:, 2);
uh1 = uh(:, 1);
uh2 = uh(:, 2);

syms B1 B2;
z1 = exp(-B1*ul1(1) - B2*ul2(1)) * sampleSL(1)*dLE;
z2 = exp(-B1*uh1(1) - B2*uh2(1)) * sampleSH(1)*dHE;
for k = 2 : lsamplelen
    z1 = z1 + exp(-B1*ul1(k) - B2*ul2(k)) * sampleSL(k)*dLE;
end
for k = 2 : hsamplelen
    z2 = z2 + exp(-B1*uh1(k) - B2*uh2(k)) * sampleSH(k)*dHE;
end
gz1 = log(z1) - log(SLintergration);
gz2 = log(z2) - log(SHintergration);
%%
b1 = zeros(1000, 1000);
b2 = zeros(1000, 1000);

array_i = 0 : 0.01 : 9.99;
array_j = 0 : 0.01 : 9.99;

for i = 1 : 1000
    i_value = array_i(i);
    parfor j = 1 : 1000
        %����˫�ֽܷ�ķ�����
        gz = (gz1 - i_value)^2 + (gz2 - array_j(j))^2;
        g = matlabFunction(gz, 'vars',{[B1, B2]});
        x = [0 0];
        options = optimset('Algorithm','Levenberg-Marquardt');
        [r, f] = lsqnonlin(g, x, [], [], options);
        b1(i, j) = r(1);
        b2(i, j) = r(2);
        disp(i);
        disp(j);
    end
end

