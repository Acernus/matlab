frames=360;
sod=1800;%源到旋转中心
sdd=2000;%源到探测器
dx=400;
dy=400;
odd=sdd-sod;
%椭球方程x^2/a^2+y^2/b^2+z^2/c^2=1
alpha=15;%倾斜角
% Rx=[1,0,0,0;0,cos(-alpha),sin(-alpha),0;0,-sin(-alpha),cos(-alpha),0;0,0,0,1];
prj=zeros(frames,400,400,'single');
%%大圆
1
for theta=0:frames-1;%旋转角度
for j=1:dy
    for i=1:dx
        alpha1=alpha*pi/180;
        theta1=theta*pi/180;
%         Rz=[cos(-theta1),sin(-theta1),0,0;-sin(-theta1),cos(-theta1),0,0;0,0,1,0;0,0,0,1];
%         s=[0,-sod,0,1];%旋转alpha角度后的源坐标
        s=[-sod*cos(alpha1)*sin(theta1),-sod*cos(alpha1)*cos(theta1),sod*sin(alpha1)];
%         detector=[i-detec_x/2,odd,j-detec_y/2,1]*Rx*Rz;
        d=[(i-dx/2)*cos(theta1)+sin(theta1)*(odd*cos(alpha1)+(dy/2-j)*sin(alpha1)),-(i-dx/2)*sin(theta1)+(odd*cos(alpha1)+(dy/2-j)*sin(alpha1))*cos(theta1), -odd*sin(alpha1)+(dy/2-j)*cos(alpha1)];
        mu=1;
        a=120;%x轴
        b=100;%y轴
        c=30;%z轴
        A = (d(1)-s(1))^2/a^2+(d(2)-s(2))^2/b^2+(d(3)-s(3))^2/c^2;
        B= 2*s(1)*(d(1)-s(1))/a^2+2*s(2)*(d(2)-s(2))/b^2+2*s(3)*(d(3)-s(3))/c^2;
        C = s(1)^2/a^2+s(2)^2/b^2+s(3)^2/c^2-1;
        dlta = B*B -4*A*C;
        if(dlta>0)
            t1  = (-B-sqrt(dlta))/(2*A);
            t2  = (-B+sqrt(dlta))/(2*A);
            x1 = t1*(d(1)-s(1))+s(1);
            x2 = t2*(d(1)-s(1))+s(1);
            y1 = t1*(d(2)-s(2))+s(2);
            y2 = t2*(d(2)-s(2))+s(2);
            z1=  t1*(d(3)-s(3))+s(3);
            z2 = t2*(d(3)-s(3))+s(3);
            prj(theta+1,j,i) =prj(theta+1,j,i)+mu*sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
        end       
        end 
end
    
end
%%大圆锥0（0,0，-30）
2
for theta=0:frames-1;%旋转角度
    for j=1:dy
    for i=1:dx
        alpha1=alpha*pi/180;
        theta1=theta*pi/180;
%         Rz=[cos(-theta1),sin(-theta1),0,0;-sin(-theta1),cos(-theta1),0,0;0,0,1,0;0,0,0,1];
%         s=[0,-sod,0,1];%旋转alpha角度后的源坐标
        s=[-sod*cos(alpha1)*sin(theta1),-sod*cos(alpha1)*cos(theta1),sod*sin(alpha1)];
%         detector=[i-detec_x/2,odd,j-detec_y/2,1]*Rx*Rz;
        d=[(i-dx/2)*cos(theta1)+sin(theta1)*(odd*cos(alpha1)+(dy/2-j)*sin(alpha1)),-(i-dx/2)*sin(theta1)+(odd*cos(alpha1)+(dy/2-j)*sin(alpha1))*cos(theta1), -odd*sin(alpha1)+(dy/2-j)*cos(alpha1)];
     r=80;%圆锥体底半径
        h=50;%圆锥体高
        mu=-1;
        x0=0;y0=0;z0=-30;
        A=(d(1)-s(1)-x0)^2+(d(2)-s(2)-y0)^2-(d(3)-s(3)-z0)^2*r*r/(h*h);
        B=2*s(1)*(d(1)-s(1)-x0)+2*s(2)*(d(2)-s(2)-y0)-2*s(3)*(d(3)-s(3)-z0)*r*r/(h*h);
        C=s(1)^2+s(2)^2-r*r*s(3)^2/(h*h);
        dlta = B*B -4*A*C;
        if(dlta>0)
            t1  = (-B-sqrt(dlta))/(2*A);
            t2  = (-B+sqrt(dlta))/(2*A);
            x1 = t1*(d(1)-s(1))+s(1);
            x2 = t2*(d(1)-s(1))+s(1);
            y1 = t1*(d(2)-s(2))+s(2);
            y2 = t2*(d(2)-s(2))+s(2);
            z1=  t1*(d(3)-s(3))+s(3);
            z2 = t2*(d(3)-s(3))+s(3);
            if(z1>z0&&z1<=z0+h&&z2>z0&&z2<=z0+h)
                prj(theta+1,j,i)=prj(theta+1,j,i)+mu*sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
                if(prj(theta+1,j,i)<0) 
                    prj(theta+1,j,i)=0;
                end
            end
        end      
    end 
    end
end 
% % 



%% 小圆锥1（-10,80，-10）
for theta=0:frames-1;%旋转角度
    for j=1:dy
    for i=1:dx
        alpha1=alpha*pi/180;
        theta1=theta*pi/180;
%         Rz=[cos(-theta1),sin(-theta1),0,0;-sin(-theta1),cos(-theta1),0,0;0,0,1,0;0,0,0,1];
%         s=[0,-sod,0,1];%旋转alpha角度后的源坐标
        s=[-sod*cos(alpha1)*sin(theta1),-sod*cos(alpha1)*cos(theta1),sod*sin(alpha1)];
%         detector=[i-detec_x/2,odd,j-detec_y/2,1]*Rx*Rz;
        d=[(i-dx/2)*cos(theta1)+sin(theta1)*(odd*cos(alpha1)+(dy/2-j)*sin(alpha1)),-(i-dx/2)*sin(theta1)+(odd*cos(alpha1)+(dy/2-j)*sin(alpha1))*cos(theta1), -odd*sin(alpha1)+(dy/2-j)*cos(alpha1)];
     r=10;%圆锥体底半径
        h=20;%圆锥体高
        mu=20;
        x0=-10;y0=80;z0=-10;
        A=(d(1)-s(1)-x0)^2+(d(2)-s(2)-y0)^2-(d(3)-s(3)-z0)^2*r*r/(h*h);
        B=2*s(1)*(d(1)-s(1)-x0)+2*s(2)*(d(2)-s(2)-y0)-2*s(3)*(d(3)-s(3)-z0)*r*r/(h*h);
        C=s(1)^2+s(2)^2-r*r*s(3)^2/(h*h);
        dlta = B*B -4*A*C;
        if(dlta>0)
            t1  = (-B-sqrt(dlta))/(2*A);
            t2  = (-B+sqrt(dlta))/(2*A);
            x1 = t1*(d(1)-s(1))+s(1);
            x2 = t2*(d(1)-s(1))+s(1);
            y1 = t1*(d(2)-s(2))+s(2);
            y2 = t2*(d(2)-s(2))+s(2);
            z1=  t1*(d(3)-s(3))+s(3);
            z2 = t2*(d(3)-s(3))+s(3);
            if(z1>z0&&z1<=z0+h&&z2>z0&&z2<=z0+h)
                prj(theta+1,j,i)=prj(theta+1,j,i)+mu*sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
            end
        end  
    end 
    end
end 
%%小圆锥2（-10，-80，-10）
for theta=0:frames-1;%旋转角度
    for j=1:dy
    for i=1:dx
        alpha1=alpha*pi/180;
        theta1=theta*pi/180;
%         Rz=[cos(-theta1),sin(-theta1),0,0;-sin(-theta1),cos(-theta1),0,0;0,0,1,0;0,0,0,1];
%         s=[0,-sod,0,1];%旋转alpha角度后的源坐标
        s=[-sod*cos(alpha1)*sin(theta1),-sod*cos(alpha1)*cos(theta1),sod*sin(alpha1)];
%         detector=[i-detec_x/2,odd,j-detec_y/2,1]*Rx*Rz;
        d=[(i-dx/2)*cos(theta1)+sin(theta1)*(odd*cos(alpha1)+(dy/2-j)*sin(alpha1)),-(i-dx/2)*sin(theta1)+(odd*cos(alpha1)+(dy/2-j)*sin(alpha1))*cos(theta1), -odd*sin(alpha1)+(dy/2-j)*cos(alpha1)];
     r=10;%圆锥体底半径
        h=20;%圆锥体高
        mu=30;
        x0=-10;y0=-80;z0=-10;
        A=(d(1)-s(1)-x0)^2+(d(2)-s(2)-y0)^2-(d(3)-s(3)-z0)^2*r*r/(h*h);
        B=2*s(1)*(d(1)-s(1)-x0)+2*s(2)*(d(2)-s(2)-y0)-2*s(3)*(d(3)-s(3)-z0)*r*r/(h*h);
        C=s(1)^2+s(2)^2-r*r*s(3)^2/(h*h);
        dlta = B*B -4*A*C;
        if(dlta>0)
            t1  = (-B-sqrt(dlta))/(2*A);
            t2  = (-B+sqrt(dlta))/(2*A);
            x1 = t1*(d(1)-s(1))+s(1);
            x2 = t2*(d(1)-s(1))+s(1);
            y1 = t1*(d(2)-s(2))+s(2);
            y2 = t2*(d(2)-s(2))+s(2);
            z1=  t1*(d(3)-s(3))+s(3);
            z2 = t2*(d(3)-s(3))+s(3);
            if(z1>z0&&z1<=z0+h&&z2>z0&&z2<=z0+h)
                prj(theta+1,j,i)=prj(theta+1,j,i)+mu*sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
            end
        end
    end 
    end
end 
%%小圆锥3（-78,76）
for theta=0:frames-1;%旋转角度
    for j=1:dy
    for i=1:dx
        alpha1=alpha*pi/180;
        theta1=theta*pi/180;
%         Rz=[cos(-theta1),sin(-theta1),0,0;-sin(-theta1),cos(-theta1),0,0;0,0,1,0;0,0,0,1];
%         s=[0,-sod,0,1];%旋转alpha角度后的源坐标
        s=[-sod*cos(alpha1)*sin(theta1),-sod*cos(alpha1)*cos(theta1),sod*sin(alpha1)];
%         detector=[i-detec_x/2,odd,j-detec_y/2,1]*Rx*Rz;
        d=[(i-dx/2)*cos(theta1)+sin(theta1)*(odd*cos(alpha1)+(dy/2-j)*sin(alpha1)),-(i-dx/2)*sin(theta1)+(odd*cos(alpha1)+(dy/2-j)*sin(alpha1))*cos(theta1), -odd*sin(alpha1)+(dy/2-j)*cos(alpha1)];
     r=10;%圆锥体底半径
        h=20;%圆锥体高
        mu=20;
        x0=-78;y0=76;z0=-10;
        A=(d(1)-s(1)-x0)^2+(d(2)-s(2)-y0)^2-(d(3)-s(3)-z0)^2*r*r/(h*h);
        B=2*s(1)*(d(1)-s(1)-x0)+2*s(2)*(d(2)-s(2)-y0)-2*s(3)*(d(3)-s(3)-z0)*r*r/(h*h);
        C=s(1)^2+s(2)^2-r*r*s(3)^2/(h*h);
        dlta = B*B -4*A*C;
        if(dlta>0)
            t1  = (-B-sqrt(dlta))/(2*A);
            t2  = (-B+sqrt(dlta))/(2*A);
            x1 = t1*(d(1)-s(1))+s(1);
            x2 = t2*(d(1)-s(1))+s(1);
            y1 = t1*(d(2)-s(2))+s(2);
            y2 = t2*(d(2)-s(2))+s(2);
            z1=  t1*(d(3)-s(3))+s(3);
            z2 = t2*(d(3)-s(3))+s(3);
            if(z1>z0&&z1<=z0+h&&z2>z0&&z2<=z0+h)
                prj(theta+1,j,i)=prj(theta+1,j,i)+mu*sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
            end
        end
    end 
    end
end 
%% 小圆锥4（78，-76，-10）
for theta=0:frames-1;%旋转角度
    for j=1:dy
    for i=1:dx
        alpha1=alpha*pi/180;
        theta1=theta*pi/180;
%         Rz=[cos(-theta1),sin(-theta1),0,0;-sin(-theta1),cos(-theta1),0,0;0,0,1,0;0,0,0,1];
%         s=[0,-sod,0,1];%旋转alpha角度后的源坐标
        s=[-sod*cos(alpha1)*sin(theta1),-sod*cos(alpha1)*cos(theta1),sod*sin(alpha1)];
%         detector=[i-detec_x/2,odd,j-detec_y/2,1]*Rx*Rz;
        d=[(i-dx/2)*cos(theta1)+sin(theta1)*(odd*cos(alpha1)+(dy/2-j)*sin(alpha1)),-(i-dx/2)*sin(theta1)+(odd*cos(alpha1)+(dy/2-j)*sin(alpha1))*cos(theta1), -odd*sin(alpha1)+(dy/2-j)*cos(alpha1)];
     r=10;%圆锥体底半径
        h=20;%圆锥体高
        mu=40;
        x0=78;y0=-76;z0=-10;
        A=(d(1)-s(1)-x0)^2+(d(2)-s(2)-y0)^2-(d(3)-s(3)-z0)^2*r*r/(h*h);
        B=2*s(1)*(d(1)-s(1)-x0)+2*s(2)*(d(2)-s(2)-y0)-2*s(3)*(d(3)-s(3)-z0)*r*r/(h*h);
        C=s(1)^2+s(2)^2-r*r*s(3)^2/(h*h);
        dlta = B*B -4*A*C;
        if(dlta>0)
            t1  = (-B-sqrt(dlta))/(2*A);
            t2  = (-B+sqrt(dlta))/(2*A);
            x1 = t1*(d(1)-s(1))+s(1);
            x2 = t2*(d(1)-s(1))+s(1);
            y1 = t1*(d(2)-s(2))+s(2);
            y2 = t2*(d(2)-s(2))+s(2);
            z1=  t1*(d(3)-s(3))+s(3);
            z2 = t2*(d(3)-s(3))+s(3);
            if(z1>z0&&z1<=z0+h&&z2>z0&&z2<=z0+h)
                prj(theta+1,j,i)=prj(theta+1,j,i)+mu*sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
            end
        end
    end 
    end
end 
%%
prj=permute(prj,[2 3 1]);
writebin('c1z1.bin',prj);

