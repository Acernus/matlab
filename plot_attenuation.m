clc;
clear all;
data = load('1_80kv_Al.txt');
x = data(:, 1);
y = data(:, 2);
c = polyfit(x, y, 1);  %进行拟合，c为2次拟合后的系数
d = polyval(c, 1.1413);  %拟合后，每一个横坐标对应的值即为d
plot(x, d, 'r');       %拟合后的曲线

plot(x, y);       %将每个点 用*画出来
hold on;
