function [res]=statistic_reconstruction(type, t, PicN)
%%
image = zeros(PicN, PicN); %ÖØ½¨Í¼Ïñ

for i = 1 : PicN
    for j = 1 : PicN
        image(i, j) = 1;
    end
end

if strcmp(type, 'ml_em')
    res = caculate_ml_em(image, t);
else if strcmp(type, 'osl_em')
    res = caculate_osl_em(image, t); 
    end
end
