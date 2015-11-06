clc;
clear all;
%%
%�������ļ��ж����ߵ����״���SL��SH��
SL = load('SL.txt');
SH = load('SH.txt');
SL = SL(:, 2)';
SH = SH(:, 2)';

[lht, lwd] = size(SL);
[hht, hwd] = size(SH);
%%
%������ҵ��ĸߵ������£��������ʵ�˥��ϵ��, size == ����������
ul = [0 1; 2 3];
uh = [4 5; 6 7];
%�ֶβ���
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

%ȡ�������ߵ��ܵ�ͶӰ����
lowprj = readbin('low.bin');
highprj = readbin('high.bin');

[height, width] = size(lowprj);
%%
%�ȶ�����������ͶӰͼ�ı���
decomposition_prj1 = zeros(height, width);
decomposition_prj2 = zeros(height, width);

for i = 1 : height
    for j = 1 : width
        %����˫�ֽܷ�ķ�����
        syms B1 B2;
        z1 = exp(-B1*ul(1, 1) - B2*ul(1, 2)) * sampleSL(1)*dE;
        z2 = exp(-B1*uh(1, 1) - B2*uh(1, 2)) * sampleSH(1)*dE;
        for k = 2 : lsamplelen
            z1 = z1 + exp(-B1*ul(i, 1) - B2*ul(i, 2)) * sampleSL(i)*dE;
        end
        for k = 2 : hsamplelen
            z2 = z2 + exp(-B1*uh(i, 1) - B2*uh(i, 2)) * sampleSH(i)*dE;
        end
        
        gz = (log(z1) - log(SLintergration) - lowprj(i, j))^2 + (log(z2) - log(SHintergration) - highprj(i, j))^2;
        g = matlabFunction(gz, 'vars',{[B1, B2]});
        options = gaoptimset('Generations', 200);
        [r, f] = ga(g, 2, [], [], [], [], 0, [], [], options);
        while f > 10^-10
            [r, f] = ga(g, 2, [], [], [], [], 0, [], [], options);
        end
        decomposition_prj1(i, j) = r(1);
        decomposition_prj2(i, j) = r(2);
    end
end

imtool(decomposition_prj1, []);
imtool(decomposition_prj2, []);


