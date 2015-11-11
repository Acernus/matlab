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

%formula = [];
for i = 1 : height
   % formularow = [];
    for j = 1 : width
        %生成双能分解的方程组
        syms B1 B2;
        z1 = exp(-B1*ul(1, 1) - B2*ul(1, 2)) * sampleSL(1)*dLE;
        z2 = exp(-B1*uh(1, 1) - B2*uh(1, 2)) * sampleSH(1)*dHE;
        for k = 2 : lsamplelen
            z1 = z1 + exp(-B1*ul(k, 1) - B2*ul(k, 2)) * sampleSL(k)*dLE;
        end
        for k = 2 : hsamplelen
            z2 = z2 + exp(-B1*uh(k, 1) - B2*uh(k, 2)) * sampleSH(k)*dHE;
        end
        
        gz = (log(z1) - log(SLintergration) - lowprj(i, j))^2 + (log(z2) - log(SHintergration) - highprj(i, j))^2;
        g = matlabFunction(gz, 'vars',{[B1, B2]});
        options = gaoptimset('Generations', 200);
        [r, f] = ga(g, 2, [], [], [], [], 0, [], [], options);
        while f > 10^-6
            [r, f] = ga(g, 2, [], [], [], [], 0, [], [], options);
        end
        decomposition_prj1(i, j) = r(1);
        decomposition_prj2(i, j) = r(2);
        %formularow = [formularow gz];
    end
%     formula = [formula; formularow];
end
%%
% [fheight, fwidth] =size(formula);
% for i = 1 : fheight
%     for j = 1 : fwidth
%         
%     end
% end


imtool(decomposition_prj1, []);
imtool(decomposition_prj2, []);


