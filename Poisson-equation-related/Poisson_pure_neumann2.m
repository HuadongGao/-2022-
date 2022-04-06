clear;clc;

for kkk = 0:4

    h  =            0.1/2^kkk;
    xh =             (0:h:1)';
    N  =           length(xh);
    
    A  = gallery('tridiag',N);

    A(1,  2) = -2;
    A(N,N-1) = -2;
    
    F  = cos(2*pi*xh);

    A        = A/h/h;
    A(1,1:2) = [1,0]; % 先把左端点定下来
    F(1    ) =     0; % u(0)定为0

    uh       = A\F  ; % 求解线性方程组
    mean_uh  = sum(uh(2:(end-1))*h + (uh(1)+uh(end))*0.5*h);
    uh       = uh - mean_uh;  % 解出来之后再拉回平均为零
    
    err = uh(1:N)- cos(2*pi*xh)/(4*pi*pi);

    err = [sqrt(sum(err(2:(end-1)).^2)*h + (err(1)^2 + err(end)^2)*0.5*h), ...
           max(abs(err))]; %两个误差l2, l-inf
    
    disp('   l2 误差                     l-inf 误差 : ')
    format longe        
    disp(err)

end
	    

plot(xh,uh(1:N),'*-b')
hold on
plot(xh,cos(2*pi*xh)/(4*pi*pi),'o-r')
legend("数值解","精确解")