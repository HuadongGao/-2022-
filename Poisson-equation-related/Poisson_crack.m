%*********************************************************************%
%| This demo is to use the 5 points central finite difference scheme  %
%| for the Poisson equation with cracks. 注意收敛阶！！！                %
%| 本程序用5点中心差分格式求解 crack boundary条件下的 Poisson 方程          %
%| - \Delta u = f(x,y),      for x \in \Omega=[-1,1]x[0,1]            %
%| \partial u/partial n = 0, for x on right half of bottom boundary   %
%| u = 0 on other parts of the boundary. f(x,y)=1 is used.            %
%|********************************************************************%
%| 该程序可以以及欢迎随意使用，更改，借鉴，分发，无需考虑版权等。              %
%| 如果你发现错误，或者你有更好的方法请联系 高华东 huadong@hust.edu.cn       %
%*********************************************************************%
for kkk = 1:7

    h  =               1/4/2^kkk; % 网格size h,逐次加密check convergence rates
    xh =               -1:h:1; % x-方向节点
    yh =                0:h:1; % y-方向节点
    Nx =           length(xh); % x-方向节点个数
    Ny =           length(yh); % y-方向节点个数
    A  = kron(speye(Ny),gallery('tridiag',Nx)) + kron(gallery('tridiag',Ny),speye(Nx)); % 生成矩阵
    A  =              1/h/h*A; % 给矩阵乘上系数
    
    itot    =                 Nx*Ny; % 总的未知的节点个数
    ibd_C   = [1,Nx,itot-Nx+1,itot]; % 总的未知的节点个数：4个角点处
    % ibd 表示 Dirichlet boundary nodes(非角点) 在所有node序列中的 index：下->左->上->右
    ibd_D   = [2:((Nx+1)/2) (Nx+1):Nx:(itot-2*Nx+1) (itot-Nx+2):(itot-1) (2*Nx):Nx:(itot-Nx)];
    ibd_all = [2:(Nx-1)     (Nx+1):Nx:(itot-2*Nx+1) (itot-Nx+2):(itot-1) (2*Nx):Nx:(itot-Nx)];
    ibd_N   = [((Nx+1)/2+1):(Nx-1)];
    
    % 首先生成 right-hand-side
    [xxh,yyh]  =    meshgrid(xh,yh); % 生成右端项, position 
    F          =     f_rhs(xxh,yyh); % 生成右端项, 取值
    F          = reshape(F',itot,1); % 生成右端项, 变成列向量
    
    % 修正 bottom 边上的矩阵系数 Neumann B.C. 只需对横向节点系数 x 0.5
    for k = 1:length(ibd_N)
        A(ibd_N(k),(ibd_N(k)-1):(ibd_N(k)+1)) = A(ibd_N(k),(ibd_N(k)-1):(ibd_N(k)+1))/2;
    end
    
    % 修正 corner 处的矩阵系数 plus 右端项 Dirichlet B.C.
    A(ibd_C,:) = A(ibd_C,:)*0;
    A(:,ibd_C) = A(:,ibd_C)*0;
    for k=1:4
        A(ibd_C(k),ibd_C(k)) = 1;
        F(ibd_C(k)         ) = 0;
    end
    
    % 修正 4 个边上的矩阵系数 plus 右端项 Dirichlet B.C.
    A(ibd_D,:) = A(ibd_D,:)*0;
    A(:,ibd_D) = A(:,ibd_D)*0;
    for k = 1:length(ibd_D)
        A(ibd_D(k),ibd_D(k)) = 1;
        F(ibd_D(k)         ) = 0;
    end
    uh  = A\F; % 计算linear system,似乎Matlab自己的\做了优化，速度极快
    
    % 生成所有 node 积分 权重:边上的节点1/2,角点节点1/4
    b_vec          = ones(Nx*Ny,1);
    b_vec(ibd_C)   =           1/4;
    b_vec(ibd_all) =           1/2;
    
    load('uh2048.mat')
    ntab = h*2048;
    ue   = uh2048(1:ntab:end,1:ntab:end); % 精确解来自于一个 5096x2048 网格解
    ue   =          reshape(ue', itot,1); % 精确解来自于一个 5096x2048 网格解
    errl2   =        sqrt(h*h*b_vec'*(uh - ue).^2); % 计算 l2-error
    errlinf =                     norm(ue-uh, inf); % 计算 l-inf error
    
    disp('   l2 误差:                      l-inf 误差 : ')
    format longe
    disp([errl2, errlinf])

end


uh = reshape(uh(1:itot),Nx,Ny)';
surf(xxh,yyh,uh)
axis([-1.1 1.1 -0.1 1.1 -0.1 0.2])

title('数值解')

% subplot(1,3,2)
% surf(xxh, yyh, f_exact(xxh,yyh) )
% title('精确解')
% 
% subplot(1,3,3)
% surf(xxh, yyh, f_exact(xxh,yyh) - uh )
% title('误差')

% 定义 right-hand-side 右端函数
function y = f_rhs(x1,x2)
    y = 1 + 0*x1;
end
