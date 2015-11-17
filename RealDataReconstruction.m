clear all;
clc;

input='E:/实验数据/filt/prj-';
len = 360;
projection = [];
for i = 1 : len
    file = readbin([input num2str(i) '.bin']);
    tmp = [projection; file(2,:)];
    projection = tmp;
end
%%
image = zeros(600, 600);
for i = 1 : 600
    for j = 1 : 600
        image(i, j) = 1;
    end
end
res = caculate_ml_em(image, projection);
