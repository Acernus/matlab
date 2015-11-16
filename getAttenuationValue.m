function [value] = getAttenuationValue(data, x)
[height, ~] = size(data);
value = data(1, 2);
for i = 1 : height
    if data(i, 1) < x
        value = data(i, 2);
    end
end