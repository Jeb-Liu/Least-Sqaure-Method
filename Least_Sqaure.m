clear;clc;
%% 要拟合的函数为周期函数
%% 符号变量
syms a b c d e;
%读数据点
file = readmatrix('data.txt');
t=file(:,1);%时间轴t
x=file(:,2);%纵轴x
f1=sym('f1',[length(x),1]);
f2=sym('f2',[length(x),1]);
f3=sym('f3',[length(x),1]);
f4=sym('f4',[length(x),1]);
f5=sym('f5',[length(x),1]);

%% 初始值
Z0 = [14.9721;0.0313879;7.85072;0.477935;0.397420];

%% 对于每个 数据点 的 各个偏导 的符号运算  式(18)
for ii=1:length(x)
	f1(ii) = (a*exp(-b*c*t(ii))*sin(sqrt(1-b^2)*c*t(ii)+d)+e-x(ii))...
              *exp(-b*c*t(ii))*sin(sqrt(1-b^2)*c*t(ii)+d);
    f2(ii) = (a*exp(-b*c*t(ii))*sin(sqrt(1-b^2)*c*t(ii)+d)+e-x(ii))...
              *( -c*t(ii)*a*exp(-b*c*t(ii))*sin(sqrt(1-b^2)*c*t(ii)+d)...
                 -b*c*t(ii)/sqrt(1-b^2)*a*exp(-b*c*t(ii))*cos(sqrt(1-b^2)*c*t(ii)+d) );
    f3(ii) = (a*exp(-b*c*t(ii))*sin(sqrt(1-b^2)*c*t(ii)+d)+e-x(ii))...
              *( -b*t(ii)*a*exp(-b*c*t(ii))*sin(sqrt(1-b^2)*c*t(ii)+d)...
                 +sqrt(1-b^2)*t(ii)*a*exp(-b*c*t(ii))*cos(sqrt(1-b^2)*c*t(ii)+d) );
    f4(ii) = (a*exp(-b*c*t(ii))*sin(sqrt(1-b^2)*c*t(ii)+d)+e-x(ii))...
              *a*exp(-b*c*t(ii))*cos(sqrt(1-b^2)*c*t(ii)+d);
    f5(ii) = (a*exp(-b*c*t(ii))*sin(sqrt(1-b^2)*c*t(ii)+d)+e-x(ii));
end

%% 符号运算求和  式(18)
F1=2*sum(f1); F2=2*sum(f2); F3=2*sum(f3); F4=2*sum(f4); F5=2*sum(f5);
%% 雅可比矩阵
Fjacobi = jacobian([F1;F2;F3;F4;F5],[a;b;c;d;e]);

%% 牛顿法求解
Zn=Z0;
for count = 1:100
    An = double( subs(Fjacobi, [a;b;c;d;e], Zn) );
    Bn = double( subs([F1;F2;F3;F4;F5], [a;b;c;d;e], Zn) );
    %% dZ
    dZ = -inv(An)*Bn;
    Zn = Zn + dZ;
    f = double( subs([F1;F2;F3;F4;F5], [a;b;c;d;e], Zn) );
    %% 到达精度停止
    if norm(dZ)<0.00001 || norm(f)<0.00001
        break;
    end
end

figure(1);
plotX=0:0.01:20;
plotY=Zn(1).*exp(-Zn(2).*Zn(3).*plotX).*sin(sqrt(1-Zn(2)).*Zn(3).*plotX+Zn(4))+Zn(5);
plot(plotX,plotY);