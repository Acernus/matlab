function PSNR = getPSNR(X, Y)
% �����ֵ�����PSNR�����������MSE
% �������YΪ�գ�����ΪX���䱾��������PSNR��MSE
if nargin<2
    D = X;
else
    if any(size(X)~=size(Y))
        error('The input size is not equal to each other!');
    end
    D = X-Y;
end
MSE = sum(D(:).*D(:))/numel(X);
fprintf('%f\n', MSE);
PSNR = 10*log10(255^2/MSE);
