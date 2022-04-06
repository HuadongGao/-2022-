clear;
format longE

for kkk = 0:7

    h  =            0.1/2^kkk;
    xh =             (0:h:1)';
    N  =           length(xh);
    
    A  = gallery('tridiag',N);

    A(1,    1:2) = [1,-1]; % 这样做可以把
    A(N,(N-1):N) = [-1,1]; % 系数矩阵对称化

    coeff_integ = [0.5, ones(1,N-2), 0.5];
    M           = [1/h/h*A,coeff_integ';coeff_integ  0];
    F           = [cos(2*pi*xh);0];
    F(1)        = F(1)/2;
    F(N)        = F(N)/2;    
    
    uh          = M\F;

    err = uh(1:N)- cos(2*pi*xh)/(4*pi*pi); % 误差函数
    err = [sqrt(sum(err(2:(end-1)).^2)*h + (err(1)^2 + err(end)^2)*0.5*h), ...
           max(abs(err))]; %两个误差l2, l-inf
    format longe
    disp('   l2 误差                     l-inf 误差 : ')
    disp(err)
end
	    

plot(xh,uh(1:N),'*-b')
hold on
plot(xh,cos(2*pi*xh)/(4*pi*pi),'o-r')
legend("数值解","精确解")