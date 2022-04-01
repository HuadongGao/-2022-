% *********************************************************************
% This demo is to use the shooting method to solve the two points
% boundary value problems: 打靶法求解两点边值问题
% y'' + 2/x*y'-2/x^2 * y = sin(ln(x))/x^2, 1 < x < 2
% with boundary values y(1) = 1, y(2) = 2.
% Here, 精确解 y(x) = c1*x + c2/x^2 - 0.3*sin(ln(x))-0.1*cos(ln(x))
% with c2 = 1/70*(8-12*sin(ln(2))-4*cos(ln(2))), c1 = 1.1-c2;
% from 北京大学出版社 《数值分析》by 张平文 李铁军； page 207 ：例题 3
% *********************************************************************


t0       = 1; % initial guess  猜测的t0,十分关键
err_flag = 1; % 记录error of 右端边界点 ---> 旗标

dt  =    0.1; % 步长
T   =      2; % 终时
xh  = 1:dt:T; % 所有时刻

yh  = zeros(2,length(xh));
zh  = zeros(2,length(xh));

while (err_flag >= 1e-13)   %解ODEs： 一个是原微分方程y，一个是导数的逼近z

    y_init  = [1,t0]'; % 初值 for y
    z_init  = [0, 1]'; % 初值 for z    
    yh(:,1) =  y_init;
    zh(:,1) =  z_init;    
    
    tc  = 1 + dt; % 目前时刻
    itc = 2;      % 目前的指标 index of tc

    while (tc < T+1e-12)   % 3级3阶Kutta格式 ---> ODE solver

        K1_y = fun_yh(xh(itc-1),yh(1,itc-1),yh(2,itc-1));
        K1_z = fun_zh(xh(itc-1),zh(1,itc-1),zh(2,itc-1));

        temp_y = yh(:,itc-1) + 0.5*dt*K1_y;
        temp_z = zh(:,itc-1) + 0.5*dt*K1_z;        
        K2_y = fun_yh(xh(itc-1)+0.5*dt,temp_y(1),temp_y(2));
        K2_z = fun_zh(xh(itc-1)+0.5*dt,temp_z(1),temp_z(2));

        temp_y = yh(:,itc-1) - dt*K1_y + 2*dt*K2_y;
        temp_z = zh(:,itc-1) - dt*K1_z + 2*dt*K2_z;        
        K3_y = fun_yh(xh(itc-1)+dt,temp_y(1),temp_y(2));
        K3_z = fun_zh(xh(itc-1)+dt,temp_z(1),temp_z(2));
        
        yh(:,itc) = yh(:,itc-1) + dt*(K1_y + 4*K2_y + K3_y)/6;
        zh(:,itc) = zh(:,itc-1) + dt*(K1_z + 4*K2_z + K3_z)/6;

        if abs(tc - T) < 1e-13
            break    
        end
        
        tc  = tc  + dt;  % 时间更新   + dt 注意：这样更新tc（大数+小数），可能会导致巨大的舍入误差累积
        itc = itc +  1;  % index更新  + 1 注意：可以用更精准的 tc = itc*dt代替
        
    end
    
    err_flag = abs(yh(1,end)-2); % 计算误差indicator(收敛离奇的快，真解近似直线?)
    disp(["在t0等于" num2str(t0) "时，右端边界误差为 " num2str(err_flag)])
    
    if (err_flag<1e-14)
        break
    end
    t0 = t0 - (yh(1,end)-2)/zh(1,end); % 用Newton迭代更新t0
    
end


plot(xh,yh(1,:),'b-*'); % 画数值解
hold ON
plot(1:(1e-3):2,y_e(1:(1e-3):2),'r-')  % 画真实解

% 计算在该步长dt下计算出来的数值解和真实解的误差(取最大)
error = norm(yh(1,:)-y_e(xh),"inf");
disp(["error in L-inf norm is :" num2str(error)])


% 定义 yh 右端函数
function y = fun_yh(x,a,b)
    y = [b;-2./x.*b + 2./x.^2*a + sin(log(x))./x.^2];
end

% 定义 zh 右端函数
function y = fun_zh(x,a,b)
    y = [b;-2./x.*b + 2./x.^2*a];
end

% 定义真实解函数
function y = y_e(x)
    c2 = 1/70*(8-12*sin(log(2))-4*cos(log(2)));
    c1 = 1.1-c2;
    y = c1*x +c2./x.^2 - 0.3*sin(log(x))-0.1*cos(log(x));
end
