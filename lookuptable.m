clc;
clear all;
%%
%从能谱文件中读出高低能谱存入SL和SH中
SL = load('80kv_Spectrum.dat');
SH = load('140kv_Spectrum.dat');
SLX = SL(:,1)';
SHX = SH(:,1)';
SL = SL(:, 2)';
SH = SH(:, 2)';

[lht, lwd] = size(SL);
[hht, hwd] = size(SH);
%%
%输入查找到的高低能谱下，两种物质的衰减系数, size == 抽样能谱数
CUL = load('0_80kv_Carbon.txt');
CUH = load('0_140kv_Carbon.txt');
FEUL = load('0_80kv_Fe.txt');
FEUH = load('0_140kv_Fe.txt');

[culh, culw] = size(CUL);
[cuhh, cuhw] = size(CUH);


ul = [];
uh = [];
for i = 1 : lwd
    cul = 0;
    feul = 0;
    for j = 1 : culh
        if SLX(i) - CUL(j, 1) > 0
            cul = CUL(j, 2);
        else
            break;
        end
    end
    for j = 1 : culh
        if SLX(i) - FEUL(j, 1) > 0
            feul = FEUL(j, 2);
        else
            break;
        end
    end
    ul = [ul;cul feul];
end



for i = 1 : hwd
    cuh = 0;
    feuh = 0;
    for j = 1 : cuhh
        if SHX(i) - CUH(j, 1) > 0
            cuh = CUH(j, 2);
        else
            break;
        end
    end
    for j = 1 : cuhh
        if SHX(i) - FEUH(j, 1) > 0
            feuh = FEUH(j, 2);
        else
            break;
        end
    end
    
    uh = [uh;cuh feuh];
end
%%

%分段采样
lsamplelen = lwd;
hsamplelen = hwd;
sampleSL = SL;
sampleSH = SH;
dE = 1;

SLintergration = 0;
SHintergration = 0;

for i = 1 : lsamplelen
    SLintergration = SLintergration + sampleSL(i) * dE;
end
for i = 1 : hsamplelen
    SHintergration = SHintergration + sampleSH(i) * dE;
end
%%

pl = [];
ph = [];

for i = 1 : 0.1 : 10
    for j = 1 : 0.1 : 10
        %生成双能分解的方程组
        syms B1 B2;
        z1 = exp(-B1*ul(1, 1) - B2*ul(1, 2)) * sampleSL(1)*dE;
        z2 = exp(-B1*uh(1, 1) - B2*uh(1, 2)) * sampleSH(1)*dE;
        for k = 2 : lsamplelen
            z1 = z1 + exp(-B1*ul(k, 1) - B2*ul(k, 2)) * sampleSL(k)*dE;
        end
        for k = 2 : hsamplelen
            z2 = z2 + exp(-B1*uh(k, 1) - B2*uh(k, 2)) * sampleSH(k)*dE;
        end
        B1 = i;
        B2 = j;
        pl = [pl, eval(subs(-log(z1) + log(SLintergration)))];
        ph = [ph, eval(subs(-log(z2) + log(SHintergration)))];
    end
end
%%
B1 = 1 : 0.1 : 10;
B2 = 1 : 0.1 : 10;


surf(B1, B2, plr);
axis tight;
colormap(hot);
shading interp;
figure;
surf(B1, B2, phr);
axis tight;
colormap(hot);
shading interp;
