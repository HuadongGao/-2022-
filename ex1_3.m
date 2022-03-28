% Huadong Gao 高华东：教学用
% 李荣华-刘播：微分方程数值解法 第四版
% page 10 第三小题：使用Euler法和改进的Euler法
% 求解 u''=-u,with u(0)=0,u'(0)=1.步长分别取h=0.1，0.05，比较算法精度。
% 我们知道其真解是 u = sin(t)。
% To solve this equation numerically,we first rewrite the equation into
% a system of equations: d {u} = {0  1}{u}, with u(0)=0,
%                        dt{v} = {-1 0}{v}       v(0)=1,
% and then apply the schemes on it.

T= 1; % 最终时间
dt = 0.1/8; % 步长
xh = 0:dt:T; % 时间节点
uh1= zeros(2,length(xh)); % 存储Euler格式计算结果
uh2= zeros(2,length(xh)); % 存储改进Euler格式结果
uh1(:,1) = [0.0 1.0]'; % 赋予初始值
uh2(:,1) = [0.0 1.0]'; % 赋予初始值

A = [0 1;-1 0];      % 方程的常系数矩阵
B = inv(eye(2)-0.5*dt*A); % Mid Point Rule 需要的逆矩阵

% 下面是主要的计算部分，我们用 while loop 进行时间推进
tc = dt;    % 现在的时间 time at current time（现在是dt，因为0时刻已经过去了）
index_tc=2; % 现在的时间所对应的 index（为了读取、存储something from/into 数组）
while tc < T+1e-12    
    uh1(:,index_tc) = uh1(:,index_tc-1)+dt*A*uh1(:,index_tc-1); % 向前Euler格式
    uh2(:,index_tc) = B*(uh2(:,index_tc-1)+0.5*dt*A*uh2(:,index_tc-1)); % 改进Euler格式     
        
    if abs(tc - T) <1e-12 % 一旦运行时间 tc 到了 T 附近
        break             % 就终止 while 循环
    end                   % break 跳出
    
    tc = tc + dt;           % 时间更新，增加dt
    index_tc = index_tc +1; % index指标增加+1
end

%下面是后处理部分，包括
%画图，计算误差等

subplot(2,2,1)
plot(xh,uh1(1,:),'-*')
hold on
plot(0:1e-4:T,sin(0:1e-4:T),'b-')
title("向前Euler法： u")
subplot(2,2,2)
plot(xh,uh1(2,:),'-o')
hold on
plot(0:1e-4:T,cos(0:1e-4:T),'g-')
title("向前Euler法： v")

subplot(2,2,3)
plot(xh,uh2(1,:),'-*')
hold on
plot(0:1e-4:T,sin(0:1e-4:T),'b-')
title("改进的Euler法： u")
subplot(2,2,4)
plot(xh,uh2(2,:),'-o')
hold on
plot(0:1e-4:T,cos(0:1e-4:T),'g-')
title("改进的Euler法： v")

% 计算并显示误差
format longE
disp(["向前Euler法误差：", num2str([abs(sin(1)-uh1(1,end)),abs(cos(1)-uh1(2,end))])])
disp(["改进的Euler法(mid point rule)误差：", num2str([abs(sin(1)-uh2(1,end)),abs(cos(1)-uh2(2,end))])])


