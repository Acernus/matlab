clc;
clear all;
%%
%�������ļ��ж����ߵ����״���SL��SH��
SL = load('newSpectrum80kV.dat');%1-80KeV
SH = load('newSpectrum160kV.dat');%1-160Kev
points = 600;
dLE = 80/points;
dHE = 160/points;
[lht, lwd] = size(SL);
[hht, hwd] = size(SH);
%%
%�ֶβ������˴���ȡ600����
lsamplelen = points;
hsamplelen = points;
sampleSL = SL;
sampleSH = SH;
%������ҵ��ĸߵ������£��������ʵ�˥��ϵ��, size == ����������
%����������ϣ��ֱ�õ��ߵ��������ʵ�˥������
L_Al = interpolant('0_80kv_Al.txt', 5);
L_Carbon = interpolant('0_80kv_Carbon.txt', 5);
H_Al = interpolant('0_160kv_Al.txt', 5);
H_Carbon = interpolant('0_160kv_Carbon.txt', 5);
ul = [];
uh = [];
for i = 0 : lwd - 1
    ulAl = getAttenuationValue(L_Al, 1 + dLE*i);
    ulCarbon = getAttenuationValue(L_Carbon, 1 + dLE*i);
    ul = [ul;ulAl ulCarbon];
end
for i = 0 : hwd - 1
    uhAl = getAttenuationValue(H_Al, 1 + dHE*i);
    uhCarbon = getAttenuationValue(H_Carbon, 1 + dHE*i);
    uh = [uh; uhAl uhCarbon];
end
%% uh document
[rows, ~] = size(uh);
fid = fopen('E:/matlab/uh.txt','wt');
for i = 1 : rows
    fprintf(fid, '%f %f\n', uh(i, 1), uh(i, 2));
end
fclose(fid);
%ul document
[rows, ~] = size(ul);
fid = fopen('E:/matlab/ul.txt','wt');
for i = 1 : rows
    fprintf(fid, '%f %f\n', ul(i, 1), ul(i, 2));
end
fclose(fid);

