clc;
clear all;
data = load('1_80kv_Al.txt');
x = data(:, 1);
y = data(:, 2);
c = polyfit(x, y, 1);  %������ϣ�cΪ2����Ϻ��ϵ��
d = polyval(c, 1.1413);  %��Ϻ�ÿһ���������Ӧ��ֵ��Ϊd
plot(x, d, 'r');       %��Ϻ������

plot(x, y);       %��ÿ���� ��*������
hold on;
