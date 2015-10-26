function writebin(file,a)
if(ndims(a)==3)
    [height width,zbot] = size(a);
    fid = fopen(file,'wb');
    for z=1:zbot
        b=squeeze(a(:,:,z));
        b=b';
    fwrite(fid,b,'float');
    end
    k = zeros(1,123);
    k = double(k);
    fwrite(fid,k,'float');
    k1 = squeeze(min(min(a(:))));
    fwrite(fid,k1,'float');
    k2 = squeeze(max(max(a(:))));
    fwrite(fid,k2,'float');
    fwrite(fid,width,'int');
    fwrite(fid,height,'int');
    Depth = zbot;
    fwrite(fid,Depth,'int');
    fclose(fid);
elseif(ismatrix(a)==2)
    [height width] = size(a);
    fid = fopen(file,'wb');
    b=a';
    fwrite(fid,b,'float');
    
    k = zeros(1,123);
    k = double(k);
    fwrite(fid,k,'float');
    k1 = (min(min(a(:))));
    fwrite(fid,k1,'float');
    k2 = (max(max(a(:))));
    fwrite(fid,k2,'float');
    fwrite(fid,width,'int');
    fwrite(fid,height,'int');
    Depth = 1;
    fwrite(fid,Depth,'int');
    fclose(fid);
end