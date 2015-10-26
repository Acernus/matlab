clc;
clear all;
%get the two models，定义全局变量
R0 = 400; R1 = 0; M = 360; N = 600; PicN = 512;
lowModelIntensity = [0.2 0.3 0.5 0.7 0.9];
highModelIntensity = [0.15 0.21 0.35 0.49 0.63];
%%
%原始模型
lowModel = create(lowModelIntensity, PicN);
highModel = create(highModelIntensity, PicN);
%%
%获取分解后的投影
g = getProjection(lowModelIntensity, R0, R1, M, N);
g1 = getProjection(highModelIntensity, R0, R1, M, N);

% In high voltage u1 = 1, u2 = 1.2, in low voltage u1 = 1.2, u2 = 1.5
lu1 = 0.2; lu2 = 0.3;
hu1 = 0.15; hu2 = 0.21;

[row, col] = size(g);
t1 = zeros(row, col);
t2 = zeros(row, col);
for i = 1 : row 
    for j = 1 : col
        t2(i, j) = (hu1 * g(i, j) - lu1 * g1(i, j)) / (hu1 * lu2 - lu1 * hu2);
        t1(i, j) = (g1(i, j) - t2(i, j) * hu2) / hu1;
    end
end
%%
%FBP重建结果
rt1 = FbpProjection(R0, R1, M, N, t1, PicN);
rt2 = FbpProjection(R0, R1, M, N, t2, PicN);
figure, imshow(rt1, []);
figure, imshow(rt2, []);
%%
%统计重建结果ml-em
ml_st1 = statistic_reconstruction('ml_em', t1, PicN);
%ml_st2 = FbpProjection(R0, R1, M, N, t2, PicN);
figure, imshow(ml_st1, []);
%figure, imshow(ml_rt2, []);
%%

% fid=fopen('t1.txt','wt');%写入的文件，各函数后面有说明?
% [m,n]=size(t1);
% for i=1:1:m
%     for j=1:1:n
%         if j==n 
%             fprintf(fid,'%g\n',t1(i,j)); 
%         else
%             fprintf(fid,'%g ',t1(i,j));
%         end
%     end
% end
% fclose(fid);
% 
% fid=fopen('t2.txt','wt');%写入的文件，各函数后面有说明?
% [m,n]=size(t2);
% for i=1:1:m
%     for j=1:1:n
%         if j==n 
%             fprintf(fid,'%g\n',t2(i,j)); 
%         else
%             fprintf(fid,'%g ',t2(i,j));
%         end
%     end
% end
% fclose(fid);
% 
%  figure, imshow(r1,[]);
%  figure, imshow(r2,[]);
 %%


%set the estimate image (the original value)
[width, height] = size(lowModel);

estimateImage = zeros(width, height);

%set the line numbers is 100
lineNums = 100;
%caculate the equation of the all x-ray lines
%create an array to store the all line estimate value
lines = zeros(lineNums);
pixelNumsInlines = zeros(lineNums);

for i = 1 : lineNums
    for j = 1 : pixelNumsInlines(i)
        
    end
end



