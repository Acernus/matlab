function [res] = interpolant(name, count)
data = load(name);

tmp = data;
for i = 1 : count
    [height, ~] = size(tmp);
    res = [];
    for j = 1 : height - 1
        res = [res; tmp(j, 1) tmp(j, 2)];
        res = [res; (tmp(j, 1) + tmp(j + 1, 1))/2 (tmp(j, 2) + tmp(j + 1, 2)) / 2];
    end
    res = [res; tmp(height, 1) tmp(height, 2)];
    tmp = res;
end
  