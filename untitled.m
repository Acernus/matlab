% a = bzeros(100, 200);
% b = bzeros(100, 200);
% for i = 1 : 100
%     for j = 1 : 200
%         a(i, j) = 1;
%         b(i, j) = 2;
%     end
% end
a = [1 2 3 4; 5 6 7 8];
b = [3 4 5 6; 7 8 9 10];
c = plusab(a, b);
