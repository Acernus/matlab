clc;
clear all;
%%

for i = 1 : 10
    for j = 1 : 12
        get_density_zeff(sprintf('E:/lm_ml_em_img2/ml_em_img_%d.bin',i), i, sprintf('E:/lm_ml_em_img1/ml_em_img_%d.bin',j), j);
    end
end
