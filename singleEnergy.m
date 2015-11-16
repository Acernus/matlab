clc;
clear all;
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


lu1 = 5.7252E-01; lu2 = 3.7010E-01;
hu1 = 3.6952E-01; hu2 = 2.9330E-01;

for i = 1 : height 
    for j = 1 : width
        decomposition_prj2(i, j) = (hu1 * lowprj(i, j) - lu1 * highprj(i, j)) / (hu1 * lu2 - lu1 * hu2);
        decomposition_prj1(i, j) = (highprj(i, j) - decomposition_prj2(i, j) * hu2) / hu1;
    end
end

writebin('single_decompositionprj1.bin', decomposition_prj1);
writebin('single_decompositionprj2.bin', decomposition_prj2);
imtool(decomposition_prj1, []);
imtool(decomposition_prj2, []);


