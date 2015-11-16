clc;
clear all;
I0 = 65520;
PrjH = readbin('Scan_H.bin');
PrjL = readbin('Scan_L.bin');
[height, width] = size(PrjH);
for i = 1 : height
    for j = 1 : width
        PrjH(i, j) = log(I0/PrjH(i, j));
        PrjL(i, j) = log(I0/PrjL(i, j));
    end
end
imshow(PrjH,[]);
figure;
imshow(PrjL,[]);

