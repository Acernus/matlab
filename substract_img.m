prj1 = readbin('PROJ_MATCH_H2.BIN');
prj2 = readbin('PROJ_MATCH_L2.BIN');
imtool(prj1, []);
imtool(prj2, []);
%����ͶӰ����ƫ��
% for i = 1 : 3
%     prj1 = [prj1(2:400,:);prj1(1,:)];
% end

sub = prj1 - prj2;
imtool(sub);