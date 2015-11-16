clc;
clear all;
%%
%从能谱文件中读出高低能谱存入SL和SH中
SL = load('newSpectrum80kV.dat');%1-80KeV
SH = load('newSpectrum160kV.dat');%1-160Kev
points = 600;
dLE = 80/points;
dHE = 160/points;
[lht, lwd] = size(SL);
[hht, hwd] = size(SH);
%%
%分段采样，此处就取600个点
lsamplelen = points;
hsamplelen = points;
sampleSL = SL;
sampleSH = SH;
%输入查找到的高低能谱下，两种物质的衰减系数, size == 抽样能谱数
%进行曲线拟合，分别得到高低能下物质的衰减曲线
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

