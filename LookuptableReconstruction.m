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
%%

[height, width] = size(lowprj);
%先定义生成两幅投影图的变量
decomposition_prj1 = zeros(height, width);
decomposition_prj2 = zeros(height, width);

load zh.mat zh;
load zl.mat zl;

B1 = 0 : 0.01 : 9.99;
B2 = 0 : 0.01 : 9.99;
[zheight, zwidth] = size(zh);
%%
for i = 1 : height
    for j = 1 : width
        tmp = [];
        for m = 1 : zheight
            for n = 1 : zwidth
                value = zh(m, n) - highprj(i, j);
                if zh(m, n) - highprj(i, j) < 10^-5
                    tmp = [tmp; m n value];
                end
            end
        end
        
        minB1 = 0;
        minB2 = 0;
        minSum = 3000;
        [rows,~]  = size(tmp); 
        for m = 1 : rows
            sum = tmp(m, 3)^2 + (zl(int32(tmp(m, 1)), int32(tmp(m, 2))) - lowprj(i, j))^2;
            if sum < minSum
                minSum = sum;
                minB1 = B1(m);
                minB2 = B2(n);
            end
        end
        disp('i:'); disp(i);
        disp('j:'); disp(j);
        decomposition_prj1(i, j) = minB1;
        decomposition_prj2(i, j) = minB2;
    end
end

writebin('10_parts_prj1.bin', decomposition_prj1);
writebin('10_parts_prj2.bin', decomposition_prj2);
imtool(decomposition_prj1, []);
imtool(decomposition_prj2, []);


