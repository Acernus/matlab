function [g]=getProjection( intensity, R0, R1, M, N)
delta_angle=2*pi/M;
d0 = intensity(1);
d1 = intensity(2);
d2 = intensity(3);
d3 = intensity(4);
d4 = intensity(5);
%以第一个模型灰度为基准，向上累加
pht=[0       0       200   200     0    d0 
     -100    100     45    45      0    (d1-d0)
     100     100     45    45      0    (d2-d0)
    -100    -100     45    45      0    (d3-d0)
     100    -100     45    45      0    (d4-d0) ];

%%
%1、投影数据获取
g=zeros(M,N,'single');%投影数据
for i=1:M%角度采样
    rmada=i*delta_angle;
    sx=R0*cos(rmada);
    sy=R0*sin(rmada);
    for j=1:N%探测器编号，从左往右
        u=(j-N/2-0.5);
        dx=u*cos(rmada+pi/2);
        dy=u*sin(rmada+pi/2);
        
        deltax=dx-sx;
        deltay=dy-sy;
        deltaz=0;%dz-sz;                    
        for k=1:5%Nu
            P=sx-pht(k,1);
            Q=sy-pht(k,2);
            R1=0;%sz-Z0(k,1);
            px1=((deltay*sin(pht(k,5))).^2+(deltax*cos(pht(k,5))).^2+2*deltax*deltay*sin(pht(k,5))*cos(pht(k,5)))/(pht(k,3).^2);
            px2=((deltay*cos(pht(k,5))).^2+(deltax*sin(pht(k,5))).^2-2*deltax*deltay*sin(pht(k,5))*cos(pht(k,5)))/(pht(k,4).^2);
            px3=0;%deltaz.^2/(c(k,1).^2);
            px=px1+px2+px3;
            q1=(2*Q*deltay*(sin(pht(k,5))).^2+2*P*deltax*(cos(pht(k,5))).^2+2*(P*deltay+Q*deltax)*sin(pht(k,5))*cos(pht(k,5)))/(pht(k,3).^2);
            q2=(2*Q*deltay*(cos(pht(k,5))).^2+2*P*deltax*(sin(pht(k,5))).^2-2*(P*deltay+Q*deltax)*sin(pht(k,5))*cos(pht(k,5)))/(pht(k,4).^2);                   
            q3=0;%2*R*deltaz/(c(k,1).^2);
            qx=q1+q2+q3;
            rx1=((Q*sin(pht(k,5))).^2+(P*cos(pht(k,5))).^2+2*P*Q*sin(pht(k,5))*cos(pht(k,5)))/(pht(k,3).^2);
            rx2=((Q*cos(pht(k,5))).^2+(P*sin(pht(k,5))).^2-2*P*Q*sin(pht(k,5))*cos(pht(k,5)))/(pht(k,4).^2);    
            rx3= 0;%R^2/(c(k,1).^2);
            rx=rx1+rx2+rx3-1;
            if qx^2-4*px*rx>0
               t1=(-qx+sqrt(qx.^2-4*px*rx))/(2*px);
               t2=(-qx-sqrt(qx.^2-4*px*rx))/(2*px);
               jx1=(dx-sx)*t1+sx;     
               jy1=(dy-sy)*t1+sy;     
               jz1=0;%(dz-sz)*t1+sz;
               jx2=(dx-sx)*t2+sx;     
               jy2=(dy-sy)*t2+sy;     
               jz2=0;%(dz-sz)*t2+sz;
               g(i,j)=g(i,j)+pht(k,6)*sqrt((jx2-jx1).^2+(jy2-jy1).^2+(jz2-jz1).^2);
            end
        end
    end
   % i  
end
