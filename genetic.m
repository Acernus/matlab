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
decomposition_prj1 = zeros(height, width);
decomposition_prj2 = zeros(height, width);
ul1 = ul(:, 1);
ul2 = ul(:, 2);
uh1 = uh(:, 1);
uh2 = uh(:, 2);
for i = 1 : height
    for j = 1 : width
        %����˫�ֽܷ�ķ�����
        syms B1 B2;
        z1 = exp(-B1*ul1(1) - B2*ul2(1)) * sampleSL(1)*dLE;
        z2 = exp(-B1*uh1(1) - B2*uh2(1)) * sampleSH(1)*dHE;
        for k = 2 : lsamplelen
            z1 = z1 + exp(-B1*ul1(k) - B2*ul2(k)) * sampleSL(k)*dLE;
        end
        for k = 2 : hsamplelen
            z2 = z2 + exp(-B1*uh1(k) - B2*uh2(k)) * sampleSH(k)*dHE;
        end
        
        gz = (log(z1) - log(SLintergration) - lowprj(i, j))^2 + (log(z2) - log(SHintergration) - highprj(i, j))^2;
        g = matlabFunction(gz, 'vars',{[B1, B2]});
        options = gaoptimset('Generations', 200);
        [r, f] = ga(g, 2, [], [], [], [], 0, [], [], options);
        decomposition_prj1(i, j) = r(1);
        decomposition_prj2(i, j) = r(2);
        disp(i);
        disp(j);
    end
end

writebin('10_parts_prj1.bin', decomposition_prj1);
writebin('10_parts_prj2.bin', decomposition_prj2);
imtool(decomposition_prj1, []);
imtool(decomposition_prj2, []);