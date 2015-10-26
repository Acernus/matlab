R0=400;%射线源到中心距离
%重建图像大小
pi = 3.1415926535;
PicN = 512;
Pic = zeros(PicN, PicN);
M = 1;
N = 600;
iterativeTime = 1;
load FBP_t1;
load FBP_t2;
sintable = zeros(M);
costable = zeros(N);

for i = 1 : M
    sintable(i) = sin((i - 1) * pi * 2 / M);
    costable(i) = cos((i - 1) * pi * 2 / M);
end

%%

belta = 2; %贝塔值
ldelta = 2;%huber函数小delta的值
for h = 1 : iterativeTime
    for m = 1 : M
        delta = zeros(N);
        L = zeros(N);%射线的积分
        LT = zeros(PicN, PicN);%每个点的积分
        LTmutilDelta = zeros(PicN, PicN);
        w = zeros(N, PicN);%加权值
        U = zeros(N);%势函数
        %计算每个像素值上穿过射线的长度累加
        for i = 1 : PicN
            for j = 1 : PicN
                for n = 1 : N
                    x0 = R0 * costable(m);
                    y0 = R0 * sintable(m);
                    x1 = -(PicN - n) / 2 * sintable(m);
                    y1 = (PicN - n) / 2 * costable(m);
                    if abs(x1 - x0) < exp(-6)
                        if x1 <= i - PicN / 2 && x1 >= i - PicN / 2 - 1
                            LT(i, j) = LT(i, j) + 1;
                        end
                    elseif abs(y1 - y0) < exp(-6)
                        if y1 >= j - PicN / 2 - 1 && y1 <= j - PicN / 2
                            LT(i, j) = LT(i, j) + 1;    
                        end
                    else
                        ix0 = i - PicN / 2;
                        iy0 = j - PicN / 2;
                        ix1 = i - PicN / 2 - 0.5;
                        iy1 = j - PicN / 2;
                        ix2 = i - PicN / 2;
                        iy2 = j - PicN / 2 - 0.5;
                        ix3 = i - PicN / 2 - 0.5;
                        iy3 = j - PicN / 2 - 0.5;
                        k = (y1 - y0) / (x1 - x0);
                        b = y1 - k * x1;
                        Join = zeros(2, 2);
                        p = 1;
                        q = 1;
                        if ix0 * k + b <= iy0 && ix0 * k + b >= iy2
                            Join(p, q) = ix0;
                            Join(p, q + 1) = k * ix0 + b; 
                            p = p + 1;
                        end
                        if ix1 * k + b <= iy1 && ix1 * k + b >= iy3
                            Join(p, q) = ix1;
                            Join(p, q + 1) = k * ix1 + b; 
                            p = p + 1;
                        end
                        if (iy0 - b) / k >= ix0 && (iy0 - b) / k <= ix1
                            Join(p, q) = (iy0 - b) / k;
                            Join(p, q + 1) = iy0; 
                            p = p + 1;
                        end
                        if (iy2 - b) / k >= ix2 && (iy2 - b) / k <= ix3
                            Join(p, q) = (iy2 - b) / k;
                            Join(p, q + 1) = iy2; 
                            p = p + 1;
                        end

                        if p == 3
                            px1 = Join(1, 1);
                            py1 = Join(1, 2);
                            px2 = Join(2, 1);
                            py2 = Join(2, 2);

                            LT(i, j) = LT(i, j) + sqrt((px1 - px2)^2 + (py1 - py2)^2);
                        end
                    end
                end
            end
        end
        %计算每条射线的累加和delta
        for n = 1 : N
            x0 = R0 * costable(m);
            y0 = R0 * sintable(m);
            x1 = -(PicN - n) / 2 * sintable(m);
            y1 = (PicN - n) / 2 * costable(m);
            count = 1;
            if abs(x1 - x0) < exp(-6)
                if floor(x1) == 0
                    x1 = 1;
                end
                for j = 1 : PicN
                    L(n) = L(n) + Pic(floor(x1), j);
                    count = count + 1;
                    w(n, j) = 1;
                end
            elseif abs(y1 - y0) < exp(-6)
                if floor(y1) == 0
                    y1 = 1;
                end
                for j = 1 : PicN
                    L(n) = L(n) + Pic(j, floor(y1));
                    count = count + 1;
                    w(n, j) = 1;
                end
            else
                for i = 1 : PicN
                    for j = 1 : PicN
                        ix0 = i - PicN / 2;
                        iy0 = j - PicN / 2;
                        ix1 = i - PicN / 2 - 0.5;
                        iy1 = j - PicN / 2;
                        ix2 = i - PicN / 2;
                        iy2 = j - PicN / 2 - 0.5;
                        ix3 = i - PicN / 2 - 0.5;
                        iy3 = j - PicN / 2 - 0.5;
                        k = (y1 - y0) / (x1 - x0);
                        b = y1 - k * x1;
                        Join = zeros(2, 2);
                        p = 1;
                        q = 1;
                        if ix0 * k + b <= iy0 && ix0 * k + b >= iy2
                            Join(p, q) = ix0;
                            Join(p, q + 1) = k * ix0 + b; 
                            p = p + 1;
                        end
                        if ix1 * k + b <= iy1 && ix1 * k + b >= iy3
                            Join(p, q) = ix1;
                            Join(p, q + 1) = k * ix1 + b; 
                            p = p + 1;
                        end
                        if (iy0 - b) / k >= ix0 && (iy0 - b) / k <= ix1
                            Join(p, q) = (iy0 - b) / k;
                            Join(p, q + 1) = iy0; 
                            p = p + 1;
                        end
                        if (iy2 - b) / k >= ix2 && (iy2 - b) / k <= ix3
                            Join(p, q) = (iy2 - b) / k;
                            Join(p, q + 1) = iy2; 
                            p = p + 1;
                        end

                        if p == 3
                            px1 = Join(1, 1);
                            py1 = Join(1, 2);
                            px2 = Join(2, 1);
                            py2 = Join(2, 2);
                            if abs(px1 - px2) == 1 || abs(py1 - py2) == 1
                                w(n, count) = 1;
                            else
                                w(n, count) = sqrt(2);
                            end
                            count = count + 1;

                            L(n) = L(n) + sqrt((px1 - px2)^2 + (py1 - py2)^2) * Pic(i, j);
                        end
                    end
                end
            end
            delta(n) = L(n) / t1(n);
        end
        %计算delta乘以长度的和
        for i = 1 : PicN
            for j = 1 : PicN
                for n = 1 : N
                    x0 = R0 * costable(m);
                    y0 = R0 * sintable(m);
                    x1 = -(PicN - n) / 2 * sintable(m);
                    y1 = (PicN - n) / 2 * costable(m);
                    if abs(x1 - x0) < exp(-6)
                        if x1 <= i - PicN / 2 && x1 >= i - PicN / 2 - 1
                            LTmutilDelta(i, j) = LTmutilDelta(i, j) + Pic(i, j) * delta(n);
                        end
                    elseif abs(y1 - y0) < exp(-6)
                        if y1 >= j - PicN / 2 - 1 && y1 <= j - PicN / 2
                            LTmutilDelta(i, j) = LTmutilDelta(i, j) + Pic(i, j) * delta(n);    
                        end
                    else
                        ix0 = i - PicN / 2;
                        iy0 = j - PicN / 2;
                        ix1 = i - PicN / 2 - 0.5;
                        iy1 = j - PicN / 2;
                        ix2 = i - PicN / 2;
                        iy2 = j - PicN / 2 - 0.5;
                        ix3 = i - PicN / 2 - 0.5;
                        iy3 = j - PicN / 2 - 0.5;
                        k = (y1 - y0) / (x1 - x0);
                        b = y1 - k * x1;
                        Join = zeros(2, 2);
                        p = 1;
                        q = 1;
                        if ix0 * k + b <= iy0 && ix0 * k + b >= iy2
                            Join(p, q) = ix0;
                            Join(p, q + 1) = k * ix0 + b; 
                            p = p + 1;
                        end
                        if ix1 * k + b <= iy1 && ix1 * k + b >= iy3
                            Join(p, q) = ix1;
                            Join(p, q + 1) = k * ix1 + b; 
                            p = p + 1;
                        end
                        if (iy0 - b) / k >= ix0 && (iy0 - b) / k <= ix1
                            Join(p, q) = (iy0 - b) / k;
                            Join(p, q + 1) = iy0; 
                            p = p + 1;
                        end
                        if (iy2 - b) / k >= ix2 && (iy2 - b) / k <= ix3
                            Join(p, q) = (iy2 - b) / k;
                            Join(p, q + 1) = iy2; 
                            p = p + 1;
                        end

                        if p == 3
                            px1 = Join(1, 1);
                            py1 = Join(1, 2);
                            px2 = Join(2, 1);
                            py2 = Join(2, 2);

                            LTmutilDelta(i, j) = LTmutilDelta(i, j) + sqrt((px1 - px2)^2 + (py1 - py2)^2) * delta(i, j);
                        end
                    end
                end
            end
        end
        %计算U
        for n = 1 : N
            x0 = R0 * costable(m);
            y0 = R0 * sintable(m);
            x1 = -(PicN - n) / 2 * sintable(m);
            y1 = (PicN - n) / 2 * costable(m);
            if abs(x1 - x0) < exp(-6)
                if floor(x1) == 0
                    x1 = 1;
                end
                for j = 1 : PicN
                    if Pic(floor(x1), j) <= ldelta
                        U(n) = U(n) + w(n, j) * Pic(floor(x1), j);
                    else
                        U(n) = U(n) + w(n, j) * ldelta;
                    end
                end
            elseif abs(y1 - y0) < exp(-6)
                if floor(y1) == 0
                    y1 = 1;
                end
                for j = 1 : PicN
                    if Pic(j, floor(y1)) <= ldelta
                        U(n) = U(n) + w(n, j) * Pic(j, floor(y1));
                    else
                        U(n) = U(n) + w(n, j) * ldelta;
                    end
                end
            else
                count = 1;
                for i = 1 : PicN
                    for j = 1 : PicN
                        ix0 = i - PicN / 2;
                        iy0 = j - PicN / 2;
                        ix1 = i - PicN / 2 - 0.5;
                        iy1 = j - PicN / 2;
                        ix2 = i - PicN / 2;
                        iy2 = j - PicN / 2 - 0.5;
                        ix3 = i - PicN / 2 - 0.5;
                        iy3 = j - PicN / 2 - 0.5;
                        k = (y1 - y0) / (x1 - x0);
                        b = y1 - k * x1;
                        p = 1;
                        if ix0 * k + b <= iy0 && ix0 * k + b >= iy2
                            
                            p = p + 1;
                        end
                        if ix1 * k + b <= iy1 && ix1 * k + b >= iy3
                            
                            p = p + 1;
                        end
                        if (iy0 - b) / k >= ix0 && (iy0 - b) / k <= ix1
                            
                            p = p + 1;
                        end
                        if (iy2 - b) / k >= ix2 && (iy2 - b) / k <= ix3
                            
                            p = p + 1;
                        end

                        if p == 3
                            if Pic(i, j) <= ldelta
                                U(n) = U(n) + w(n, count) * Pic(i, j);
                            else
                                U(n) = U(n) + w(n, count) * ldelta;
                            end
                            count = count + 1;
                        end
                    end
                end
            end
        end
        %修正像素值
        for n = 1 : N
            x0 = R0 * costable(m);
            y0 = R0 * sintable(m);
            x1 = -(PicN - n) / 2 * sintable(m);
            y1 = (PicN - n) / 2 * costable(m);
            if abs(x1 - x0) < exp(-6)
                if floor(x1) == 0
                    x1 = 1;
                end
                for j = 1 : PicN
                    cj = LTmutilDelta(floor(x1), j) / (LT(floor(x1), j) + belta * U(n));
                    Pic(floor(x1), j) = Pic(floor(x1), j) * cj;
                end
            elseif abs(y1 - y0) < exp(-6)
                if floor(y1) == 0
                    y1 = 1;
                end
                for j = 1 : PicN
                    cj = LTmutilDelta(j, floor(y1)) / (LT(j, floor(y1)) + belta * U(n));
                    Pic(j, floor(y1)) = Pic(j, floor(y1)) * cj;
                end
            else
                for i = 1 : PicN
                    for j = 1 : PicN
                        ix0 = i - PicN / 2;
                        iy0 = j - PicN / 2;
                        ix1 = i - PicN / 2 - 0.5;
                        iy1 = j - PicN / 2;
                        ix2 = i - PicN / 2;
                        iy2 = j - PicN / 2 - 0.5;
                        ix3 = i - PicN / 2 - 0.5;
                        iy3 = j - PicN / 2 - 0.5;
                        k = (y1 - y0) / (x1 - x0);
                        b = y1 - k * x1;
                        p = 1;
                        if ix0 * k + b <= iy0 && ix0 * k + b >= iy2
                            
                            p = p + 1;
                        end
                        if ix1 * k + b <= iy1 && ix1 * k + b >= iy3
                            
                            p = p + 1;
                        end
                        if (iy0 - b) / k >= ix0 && (iy0 - b) / k <= ix1
                            
                            p = p + 1;
                        end
                        if (iy2 - b) / k >= ix2 && (iy2 - b) / k <= ix3
                            
                            p = p + 1;
                        end

                        if p == 3
                            cj = LTmutilDelta(i, j) / (LT(i, j) + belta * U(n));
                            Pic(i, j) = Pic(i, j) * cj;
                        end
                    end
                end
            end
        end
    end
end

%%
figure, imshow(Pic,[]);
