clear;
dbstop if error;

%% Load Image:
img = imread('~/Dropbox/lasso/nov21/horline.png');

[L, C, D, N] = init_4way_sparse(img);

%d  = zeros([3*N, 1]); d(3:3:end) = -1; % d'*x <= 0 <=> c_i >= 0 
d  = zeros([N, 3*N]);
for i=1:N
   d(i, 3*i) = -1; 
end

xi = 10.0;  % noramally called lambda
mu = 1.0;

%% The problem is:
%% min (0.5*z'(A'/Q*A)z+(b'/Q*A+l)z  
%% where z = [u_{1-4} v w];
%% x_op = -Q\(Az+b)

%%  y
y = zeros([3*N, 1]);
y(1:3:end) = -1; y(3:3:end) =0.5;

%% these are independent of y
Q = L+mu*D; % Q = sparse(Q); 
Q_inv = 0.5*(eye(size(Q))/Q+ Q\eye(size(Q)));
l = zeros([1, 14*N]); l(end-N:end)=1;
Aeq = []; beq = zeros(4*2*N,1);
for i=1:4*N
   Aeq = blkdiag(Aeq, [1 -1 0; 0 1 -1]); 
end
Aeq = [Aeq zeros(2*4*N, 2*N)];

%% Make Y matrix
Y = zeros([3*N, N]);
for i=1:N
   Y(3*(i-1)+1: 3*(i-1)+2, i) = y(3*(i-1)+1: 3*(i-1)+2);
end
b = -mu*D*y;

A = [C(:,:,1)', C(:,:,2)', C(:,:,3)', C(:,:,4)', d', Y];
A = sparse(A);

%% let H=(A'/Q*A), f=(b'/Q*A+l)
H = A'*Q_inv*A;
f = b'*Q_inv*A+l;

lb = zeros([14*N,1]); ub = lb;
lb(1:12*N) = -xi; lb(13*N+1:end) = [];
ub(1:12*N) =  xi; ub(12*N+1:end) = [];

% x = quadprog(H,f,A,b,Aeq,beq,lb,ub)
opts = optimoptions('quadprog','Algorithm','active-set','Display','off');
[z_op, fval] = quadprog(H, f, [], [], Aeq, beq, lb, ub, zeros([14*N, 1]), opts);
x_op = -Q_inv*(A*z_op+b);
