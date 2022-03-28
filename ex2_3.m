% Huadong Gao 高华东：教学用
% 李荣华-刘播：微分方程数值解法 第四版
% page 22 第三小题：使用3阶Adams-Bashforth（显）以及Adams-Moulton（隐）方法
% 求解 u'=-5u,u(0)=1.步长分别取h=0.1，0.05，比较算法精度。
% 首先，我们知道其真解是 u = exp(-5t)。At t=1, u= exp(-5).
% AB3: u_(n+1) = u_n + h/12*(23*f_n   -16*f_(n-1)+5*f_(n-2))
% AM3: u_(n+1) = u_n + h/12*( 5*f_(n+1)+8*f_(n)    -f_(n-1))
%

T= 1; % 最终时间
dt = 0.1; % 步长
xh = 0:dt:T; % 时间节点
uh1= zeros(length(xh),1); % 存储AB3格式计算结果
uh2= zeros(length(xh),1); % 存储AM3格式计算结果

%%% 对初始的几层数值解赋予初值：多步方法对前几步的数值解都有依赖AB3依赖前面3层的解，AM3依赖2层。
uh1(1) = 1.0;             % 赋予初始值
uh1(2) = exp(-5*dt);      % 赋予初始值
uh1(3) = exp(-5*2*dt);    % 赋予初始值
uh1(2) = (uh1(1)+0.5*dt*(-5*uh1(1)))/(1+0.5*dt*5); % 梯形公式凑合--如何提供足够精度的初始数值解？
uh1(3) = (uh1(2)+dt/12*(8*(-5*uh1(2))-(-5*uh1(1))))/(1-dt/12*5*(-5)); % 使用隐式AM3算第3步

uh2(1) = 1.0;             % 赋予初始值
uh2(2) = exp(-5*dt);      % 赋予初始值:作弊--我们直接用真实解代替
uh2(2) = (uh2(1)+0.5*dt*(-5*uh2(1)))/(1+0.5*dt*5); % 梯形公式凑合--如何提供足够精度的初始数值解？
uh2(3) = (uh2(2)+dt/12*(8*(-5*uh2(2))-(-5*uh2(1))))/(1-dt/12*5*(-5)); % 使用AM3算第3步，一起循环


% 下面是主要的计算部分，我们用while loop进行时间推进
tc = 3*dt;    % 现在的时间 time at current time（现在是3*dt，因为好几个步长已经过去了）
itc=4;        % 现在的时间所对应的 index（为了读取、存储something from/into 数组）
while tc < T+1e-12    

    uh1(itc)= uh1(itc-1)+dt/12*(23*(-5*uh1(itc-1))-16*(-5*uh1(itc-2)) +5*(-5*uh1(itc-3))); % AB3
    uh2(itc)= (uh2(itc-1)+dt/12*(8*(-5*uh2(itc-1))-(-5*uh2(itc-2))))/(1-dt/12*5*(-5)); % AM3  
        
    if abs(tc - T) <1e-12 % 一旦运行时间 tc 到了 T 附近
        break             % 就终止 while 循环
    end                   % break 跳出
    
    tc = tc + dt;         % 时间更新，增加dt
    itc = itc +1;         % index指标增加+1
end

%下面是后处理部分，包括
%画图，计算误差等

subplot(1,2,1)
plot(xh,uh1,'-*')
hold on
plot(0:1e-4:T,exp(-5*(0:1e-4:T)),'r-')
title("3阶Adams-Bashforth格式")

subplot(1,2,2)
plot(xh,uh2,'-o')
hold on
plot(0:1e-4:T,exp(-5*(0:1e-4:T)),'r-')
title("3阶Adams-Moulton格式")

% 计算并显示误差
format longE
disp([abs(exp(-5)-uh1(end)),abs(exp(-5)-uh2(end))])


