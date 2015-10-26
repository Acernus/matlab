function [Model]=create(intensity, size)
d0 = intensity(1);
d1 = intensity(2);
d2 = intensity(3);
d3 = intensity(4);
d4 = intensity(5);
Model=zeros(size, size);
width = size;
height = size;
mid = width / 2;
rmax = mid * mid * 0.49;
rmin = mid * mid * 0.0225;
for i = 1 : width
    for j = 1 : height
        if((i - mid) * (i - mid) + (j - mid) * (j - mid) > rmax)
            Model(i, j) = 0;
        else
            Model(i, j) = d0;
        end
    end
end

for i = 1 : width
    for j = 1 : height
        if((i - 0.7 * mid) * (i - 0.7 * mid) + (j - 0.7 * mid) * (j - 0.7 * mid) <= rmin)
            Model(i, j) = d1;
        end
    end
end

for i = 1 : width
    for j = 1 : height
        if((i - 0.7 * mid) * (i - 0.7 * mid) + (j - 1.3* mid) * (j - 1.3 * mid) <= rmin)
            Model(i, j) = d2;
        end
    end
end


for i = 1 : width
    for j = 1 : height
        if((i - 1.3 * mid) * (i - 1.3 * mid) + (j - 0.7 * mid) * (j - 0.7 * mid) <= rmin)
            Model(i, j) = d3;
        end
    end
end

for i = 1 : width
    for j = 1 : height
        if((i - 1.3 * mid) * (i - 1.3 * mid) + (j - 1.3 * mid) * (j - 1.3 * mid) <= rmin)
            Model(i, j) = d4;
        end
    end
end