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
%%
sumH = 0;
sumL = 0;
single_uh_carbon = 0;
single_ul_carbon = 0;
single_uh_al = 0;
single_ul_al = 0;

[~,len] = size(sampleSH); 
for i = 1 : len
    sumH = sumH + sampleSH(i);
    sumL = sumL + sampleSL(i);
end

for i = 1 : len
    single_uh_carbon = single_uh_carbon + sampleSH(i) / sumH * uh(i, 1);
    single_uh_al = single_uh_al + sampleSH(i) / sumH * uh(i, 2);
    
    single_ul_carbon = single_ul_carbon + sampleSL(i) / sumL * ul(i, 1);
    single_ul_al = single_ul_al + sampleSL(i) / sumL * ul(i, 2);
end
%%
%取出两幅高低能的投影数据
highprj = readbin('PROJ_MATCH_H2.BIN');
lowprj = readbin('PROJ_MATCH_L2.BIN');
%高能投影三次偏移
for i = 1 : 3
    highprj = [highprj(2:400,:);highprj(1,:)];
end

[height, width] = size(lowprj);

%%
%先定义生成两幅投影图的变量
decomposition_prj1 = zeros(height, width);
decomposition_prj2 = zeros(height, width);


lu1 = single_ul_carbon; lu2 = single_ul_al;
hu1 = single_uh_carbon; hu2 = single_uh_al;

for i = 1 : height 
    for j = 1 : width
        decomposition_prj2(i, j) = (hu1 * lowprj(i, j) - lu1 * highprj(i, j)) / (hu1 * lu2 - lu1 * hu2);
        decomposition_prj1(i, j) = (highprj(i, j) - decomposition_prj2(i, j) * hu2) / hu1;
    end
end

imtool(decomposition_prj1, []);
imtool(decomposition_prj2, []);


