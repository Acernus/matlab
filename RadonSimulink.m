%%
clc;
clear all;
%get the two models
% size = 200 * 200
lowModel = create(0.2, 0.3, 0.5, 0.7, 0.9, 200);
highModel = create(0.15, 0.21, 0.35, 0.49, 0.63, 200);
figure, imshow(lowModel, [0 1]);
figure, imshow(highModel, [0 1]);

% get the two projections
lowProjection = radon(lowModel, 0 : 359);
figure, imshow(lowProjection, []);
highProjection = radon(highModel, 0 : 359);
figure, imshow(highProjection, []);

% In high voltage u1 = 1, u2 = 1.2, in low voltage u1 = 1.2, u2 = 1.5
lu1 = 0.2; lu2 = 0.3;
hu1 = 0.15; hu2 = 0.21;

% calculate seperated projection and reconstruction
g = lowProjection;
g1 = highProjection;
[row, col] = size(g);
t1 = zeros(row, col);
t2 = zeros(row, col);
for i = 1 : row 
    for j = 1 : col
        t2(i, j) = (hu1 * g(i, j) - lu1 * g1(i, j)) / (hu1 * lu2 - lu1 * hu2);
        t1(i, j) = (g1(i, j) - t2(i, j) * hu2) / hu1;
    end
end

r1 = iradon(t1, 0 : 359);
r2 = iradon(t2, 0 : 359);
figure, imshow(r1, [0 1]);
figure, imshow(r2, [0 1]);
%%