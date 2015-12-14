function [imgpe, imgzeff] = get_density_zeff(filename1, num1, filename2, num2)
img1 = readbin(filename1);
img2 = readbin(filename2);

%density of carbon
pe1 = 2 * 2.0 * 6 / 12;
pe2 = 2 * 2.7 * 13 /27; 
[height, width] = size(img1);


for i = 1 : height
    for j = 1 : width
        if img1(i, j) <= 0
            img1(i, j) = 0;
        end
        if img2(i, j) <= 0
            img2(i, j) = 0;
        end
    end
end

imgpe = zeros(height, width);
imgzeff = zeros(height, width);

for i = 1 : height
    for j = 1 : width
        imgpe(i, j) = img1(i, j) * pe1 + img2(i, j) * pe2; 
        imgzeff(i, j) = ((img1(i, j) * pe1 * (6^3.5) + img2(i, j) * pe2 * (13^3.5)) / imgpe(i, j))^(1/3.5);
    end
end
 writebin(sprintf('E:/compare/ml_density_%d_%d.bin',num1, num2),imgpe);
 writebin(sprintf('E:/compare/ml_zeff_%d_%d.bin',num1, num2), imgzeff);
 
 %%
 img = readbin('E:/compare/ML_ZEFF_8_9.bin');
 [a, b] = size(img);
 for i = 1 : a
    for j = 1 : b
        if img(i, j) > 12.8
            img(i, j) = 6;
        end
    end
 end
 imtool(img, []);
