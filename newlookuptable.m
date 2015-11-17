clc;
clear all;
%%
%从能谱文件中读出高低能谱存入SL和SH中
SL = load('newSpectrum80kV.dat');%1-80KeV
SH = load('newSpectrum160kV.dat');%1-160Kev
points = 10;
dLE = 80/points;
dHE = 160/points;
[lht, lwd] = size(SL);
[hht, hwd] = size(SH);
%%
%分段采样，此处就取600个点
lsamplelen = points;
hsamplelen = points;
sampleSL = [];
sampleSH = [];
for i = 0 : points -1 
    sampleSL = [sampleSL SL(1 + i*dLE)];
    sampleSH = [sampleSH SH(1 + i*dHE)];
end
%输入查找到的高低能谱下，两种物质的衰减系数, size == 抽样能谱数
%进行曲线拟合，分别得到高低能下物质的衰减曲线
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
%先定义生成两幅投影图的变量
ul1 = ul(:, 1);
ul2 = ul(:, 2);
uh1 = uh(:, 1);
uh2 = uh(:, 2);
%formula = [];
zh = zeros(1000, 1000);
zl = zeros(1000, 1000);
count = 1;
for i = 0 : 0.01 : 9.99
    for j = 0 : 0.01 : 9.99
        z1 = 0;
        z2 = 0;
        for k = 1 : lsamplelen
            z1 = z1 + exp(-i*ul1(k) - j*ul2(k)) * sampleSL(k)*dLE;
            z2 = z2 + exp(-i*uh1(k) - j*uh2(k)) * sampleSH(k)*dHE;
        end
        zl(int32(i*100 + 1), int32(j*100 + 1)) = -log(z1) + log(SLintergration);
        zh(int32(i*100 + 1), int32(j*100 + 1)) = -log(z2) + log(SHintergration);
        disp('i:');disp(i);
        disp('j:');disp(j);
    end
end
%%
B1 = 0 : 0.01 : 9.99;
B2 = 0 : 0.01 : 9.99;
save zh.mat zh;
save zl.mat zl;
%%
surf(B1, B2, zh);
axis tight;
colormap(hot);
shading interp;
figure;
surf(B1, B2, zl);
axis tight;
colormap(hot);
shading interp;


