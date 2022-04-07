%*********************************************************************%
%| This demo is to use the 5 points central finite difference scheme  %
%| for the Poisson equation with Dirichlet boundary condition         %
%| 本程序用5点中心差分格式求解 Dirichlet boundary条件下的 Poisson 方程      %
%| - \Delta u = f(x,y),      for x \in \Omega=[0,1]x[0,1]             %
%|          u = g(x,y),      for x on boundary                        %
%| hx 和 hy 相互独立,可以自由选取,注意当中对称化矩阵的使用                   %
%| 程序中有验证误差过程,故外层有for kkk 循环                               %
%|********************************************************************%
%| 该程序可以(欢迎)随意使用，更改，借鉴，分发，无需考虑版权等。                %
%| 如果你发现错误，或者你有更好的方法请联系 高华东 huadong@hust.edu.cn       %
%*********************************************************************%

for kkk = 1:7
    
    hx =            1/2/2^kkk; % 网格size h,逐次加密 check convergence rates
    hy =                   hx; % 网格size h,逐次加密 check convergence rates    
    xh =               0:hx:1; % x-方向节点
    yh =               0:hy:1; % y-方向节点
    Nx =           length(xh); % x-方向节点个数
    Ny =           length(yh); % y-方向节点个数
    A  = kron(speye(Ny),gallery('tridiag',Nx))/hx/hx ...
         + kron(gallery('tridiag',Ny),speye(Nx))/hy/hy; % 生成矩阵
    
    itot =                Nx*Ny; % 总的未知的节点个数
    ibd_C=[1,Nx,itot-Nx+1,itot]; % 总的未知的节点个数：4个角点处
    % ibd_D 表示 Dirichlet boundary 结点(非角点) 在所有结点序列中的 index：下->左->上->右
    ibd_D=[2:(Nx-1) (Nx+1):Nx:(itot-2*Nx+1) (itot-Nx+2):(itot-1) (2*Nx):Nx:(itot-Nx)];
    
    % 首先生成 right-hand-side
    [xxh,yyh] =                             meshgrid(xh,yh); % 结点
    node      = [reshape(xxh',itot,1) reshape(yyh',itot,1)]; % 列结点
    F         =                  f_rhs(node(:,1),node(:,2)); % 生成右端项, 取值
    
    % 修正 corner 处的矩阵系数 plus 右端项 Dirichlet B.C.
    A(ibd_C,:) = A(ibd_C,:)*0;
    A(:,ibd_C) = A(:,ibd_C)*0;
    for k=1:4
        A(ibd_C(k),ibd_C(k)) = 1;
        F(ibd_C(k)         ) = f_g(node(ibd_C(k),1),node(ibd_C(k),2));
    end
    
    % 修正 4 个边上的矩阵系数 plus 右端项 Dirichlet B.C.
    A(ibd_D,:) = A(ibd_D,:)*0;
    for k = 1:length(ibd_D)
        A(ibd_D(k),ibd_D(k)) = 1;
        F(ibd_D(k)         ) = f_g(node(ibd_D(k),1),node(ibd_D(k),2));
    end
    
    ibd_Db =                2:(Nx-1); % 为了使得
    ibd_Dl = (Nx+1):Nx:(itot-2*Nx+1); % 矩阵对称
    ibd_Dt =    (itot-Nx+2):(itot-1); % 花费一些功夫
    ibd_Dr =     (2*Nx):Nx:(itot-Nx); % 是否有更好的方式?
    % 对称化处理系数矩阵,先上下边界
    for k = 1:length(ibd_Db)
        F(ibd_Db(k)+Nx) = F(ibd_Db(k)+Nx) - ...
            A(ibd_Db(k)+Nx,ibd_Db(k))*f_g(node(ibd_Db(k),1),node(ibd_Db(k),2));    
        A(ibd_Db(k)+Nx,ibd_Db(k)) = 0;
    
        F(ibd_Dt(k)-Nx) = F(ibd_Dt(k)-Nx) - ...
            A(ibd_Dt(k)-Nx,ibd_Dt(k))*f_g(node(ibd_Dt(k),1),node(ibd_Dt(k),2));    
        A(ibd_Dt(k)-Nx,ibd_Dt(k)) = 0;    
    end
    % 对称化处理系数矩阵,后左右边界
    for k = 1:length(ibd_Dl)
        F(ibd_Dl(k)+1) = F(ibd_Dl(k)+1) - ...
            A(ibd_Dl(k)+1,ibd_Dl(k))*f_g(node(ibd_Dl(k),1),node(ibd_Dl(k),2));    
        A(ibd_Dl(k)+1,ibd_Dl(k)) = 0;
    
        F(ibd_Dr(k)-1) = F(ibd_Dr(k)-1) - ...
            A(ibd_Dr(k)-1,ibd_Dr(k))*f_g(node(ibd_Dr(k),1),node(ibd_Dr(k),2));    
        A(ibd_Dr(k)-1,ibd_Dr(k)) = 0;    
    end
    
    uh  = A\F; % 计算linear system,似乎Matlab自己的\做了优化，速度极快
    
    b_vec        = ones(Nx*Ny,1); % 生成所有 node 积分 权重
    b_vec(ibd_C) =           1/4; % 边上的节点1/2
    b_vec(ibd_D) =           1/2; %   角点节点1/4

    ue      =        f_u(node(:,1),node(:,2)); % 获取真实解
    errl2   = sqrt(hx*hy*b_vec'*(uh - ue).^2); % 计算 l2-error
    errlinf =                norm(ue-uh, inf); % 计算 l-inf error
    
    disp([' l2 误差: ' '网格 nx = ' num2str(1/hx) '  ny = ' num2str(1/hy) '  l-inf 误差 : '])
    format longe
    disp([errl2, errlinf])

end

%**********************************
% 画图和后处理部分
%**********************************
uh = reshape(uh(1:itot),Nx,Ny)';
ue = reshape(ue(1:itot),Nx,Ny)';

subplot(1,3,1)
surf(xxh,yyh,uh)
title('数值解')

subplot(1,3,2)
surf(xxh,yyh,ue)
title('精确解')

subplot(1,3,3)
surf(xxh, yyh, ue - uh)
title('误差')

%**********************************
% 定义 真解, boundary data, 右端项
%**********************************
function y = f_u(x1,x2)   % 定义 exact solution 真实解
    y = exp(4*x1) + cos(4*pi*x2);
end

function y = f_g(x1,x2)   % 定义 boundary data 边界数据
    y = exp(4*x1) + cos(4*pi*x2);
end

function y = f_rhs(x1,x2) % 定义 right-hand-side 右端函数
    y = -16*exp(4*x1) + 16*pi*pi*cos(4*pi*x2);
end