%************************************************************************%
%| This demo is to use the 5 points central finite difference scheme     %
%| to solve the Poisson equation with zero Dirichlet boundary condition  %
%| on an L-shape domain                                                  %
%|                                                                       %
%|               - \Delta u=1,  x \in \Omega=[-1,1]x[-1,1]-[0,1]x[-1,0]  %
%|                        u=0 on the boundary.                           %
%|                                                                       %
%| There is no exact solution(with simple form)                          %
%|***********************************************************************%
%| This short matlab code is as demo code in my class 2022 Spring.       %
%| If you find mistakes or have a better idea in the implementation,     %
%| please send an email to Huadong GAO (Email:huadong@hust.edu.cn)       %
%************************************************************************%
clc; clear; format longE;

h  =                  1/32; % mesh size in [0,1]
xh =              (0:h:1)'; % nodes in [0,1]
N  =          length(xh)-1; % mesh number in [0,1]

Nt = N*(N-1)+(2*N-1)*(N-1); % the total unknown nodes

Ahx  = gallery('tridiag', N-1  ); % the difference matrix in x,square
Ahy  = gallery('tridiag', N    ); % the difference matrix in y,square 

Bhx  = gallery('tridiag', 2*N-1); % the difference matrix in y,rectangle
Bhy  = gallery('tridiag', N-1  ); % the difference matrix in y,rectangle

I1   = speye(N-1  ); % some identity matrix
I11  = speye(N    ); % some identity matrix
I2   = speye(2*N-1); % some identity matrix

K1 = kron(I11, Ahx) + kron(Ahy, I1); % difference matrix for square
K2 = kron(I1 , Bhx) + kron(Bhy, I2); % difference matrix for rectangle

A = sparse(Nt,Nt);  % this is the difference matrix for the whole domain
A(1:(N*(N-1)),1:(N*(N-1))) = K1;
A((N*(N-1)+1):end,(N*(N-1)+1):end) = K2;
A(((N-1)*(N-1)+1):(N*(N-1)),(N*(N-1)+1):(N*(N-1)+N-1)) = -I1;
A((N*(N-1)+1):(N*(N-1)+N-1),((N-1)*(N-1)+1):(N*(N-1))) = -I1;

uh = A\(h*h*ones(Nt,1)); % We simply use Matlab's  solver.

% plot figure: the square part
usquare = reshape(uh(1:(N*(N-1))), N-1,N)';
[xxh, yyh] = meshgrid((-1):h:0, (-1):h:0);

surf(xxh,yyh,[zeros(1,N+1); zeros(N,1) usquare zeros(N,1) ])
hold on

% plot figure: the rectangle part
usquare = reshape(uh(1:(N*(N-1))), N-1,N)';
[xxh, yyh] = meshgrid((-1):h:1, 0:h:1);
urectangle1 = [uh(((N-1)*(N-1)+1):(N*(N-1)))' zeros(1,N)];
urectangle2 = reshape(uh((N*(N-1)+1):end), 2*N-1,N-1)';
urectangle = [urectangle1;urectangle2];

[xxh, yyh] = meshgrid((-1):h:1, 0:h:1);
surf(xxh,yyh,[zeros(N,1) urectangle zeros(N,1) ; zeros(1,2*N+1)])

