% Huadong Gao 高华东：教学用
% 李荣华-刘播：微分方程数值解法 第四版
% page 10 第一小题：使用Euler法和改进的Euler法
% 求解 u'=-5u,u(0)=1.步长分别取h=0.1，0.05，比较算法精度。
% 首先，我们知道其真解是 u = exp(-5t)。At t=1, u= exp(-5).


T= 1; % 最终时间
dt = 0.1; % 步长
xh = 0:dt:T; % 时间节点
uh1= zeros(length(xh),1); % 存储Euler格式计算结果
uh2= zeros(length(xh),1); % 存储改进Euler格式结果
uh1(1) = 1.0; % 赋予初始值
uh2(1) = 1.0; % 赋予初始值

% 下面是主要的计算部分，我们用while loop进行时间推进
tc = dt;    % 现在的时间 time at current time（现在是dt，因为0时刻已经过去了）
index_tc=2; % 现在的时间所对应的 index（为了读取、存储something from/into 数组）
while tc < T+1e-12    
    uh1(index_tc)= uh1(index_tc-1)+dt*(-5*uh1(index_tc-1)); % 向前Euler格式
    uh2(index_tc)= (uh2(index_tc-1)+0.5*dt*(-5*uh2(index_tc-1)))/(1+0.5*dt*5); % 改进Euler格式     
        
    if abs(tc - T) <1e-12 % 一旦运行时间 tc 到了 T 附近
        break             % 就终止 while 循环
    end                   % break 跳出
    
    tc = tc + dt;           % 时间更新，增加dt
    index_tc = index_tc +1; % index指标增加+1
end

%下面是后处理部分，包括
%画图，计算误差等

subplot(1,2,1)
plot(xh,uh1,'-*')
hold on
plot(0:1e-4:T,exp(-5*(0:1e-4:T)),'b-')
title("向前Euler格式")

subplot(1,2,2)
plot(xh,uh2,'-o')
hold on
plot(0:1e-4:T,exp(-5*(0:1e-4:T)),'g-')
title("改进的Euler格式（Mid Point Rule）")

% 计算并显示误差
format longE
disp([abs(exp(-5)-uh1(end)),abs(exp(-5)-uh2(end))])


