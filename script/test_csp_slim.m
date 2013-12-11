dbstop if error

%% Load Image
img = imread('~/Dropbox/lasso/nov21/horline.png');

[A, B, b, D, N, ~] = init_4way_sparse(img);

%% C*x<=0 >= x_ic >=0
C = zeros([N, 3*N]);
row = 1:N; col = 3*(1:N); idx = sub2ind(size(C), row, col);
C(idx) = -1;

%% Control parameters
lambda = 1e0;
mu     = 1.0;

%% Params indepdent of y: Q and Q_inv
Q = A + mu*D;
Q_inv = eye(size(Q))/Q; Q_inv = 0.5*(Q_inv+Q_inv');

Aeq = []; beq = zeros(4*2*N,1);
for i=1:N
   Aeq = blkdiag(Aeq, [1 -1 0; 0 1 -1], [1 -1 0; 0 1 -1], [1 -1 0; 0 1 -1], [1 -1 0; 0 1 -1]);
end
Aeq = [Aeq zeros(2*4*N, 2*N)];

M = 14; % # of param blocks, total params = M*N
%% Params dependent on y:
l = zeros([1, M*N]); l(end-N:end) = 1; % l'p = \sum(w_i)
P = [B(:,:,1)', B(:,:,2)', B(:,:,3)', B(:,:,4)', C', zeros(3*N,N)];

%% Upper and lower bounds:
lb = zeros([M*N, 1]); ub = lb;
lb(1:(M-2)*N) = -lambda; lb((M-1)*N+1:end) = -Inf;
up(1:(M-2)*N) =  lambda; ub((M-2)*N+1:end) = +Inf;

%% The initial guess: y
y = zeros([3*N, 1]);
y(1:3:end) = sign(randn(N,1)).*rand(N,1);
y(2:3:end) = sign(randn(N,1)).*sqrt(1-y(1:3:end).^2);
y(3:3:end) = rand(N, 1);

%% Solving the problem:
f_prev = 1e5;
Y = zeros([3*N, N]); % place holder for y-dependent data
for it=1:100
   row = 3*(1:N); col = 1:N;
   Y(sub2ind(size(Y), row-2, col)) = y(3*(1:N)-2);
   Y(sub2ind(size(Y), row-1, col)) = y(3*(1:N)-1);
   P(:, end-N+1:end) = Y;
   q = -mu*D*y;
   
   H = P'*Q_inv*P; H = sparse(0.5*(H+H'));
   f = q'*Q_inv*P+l;
   
   % x = quadprog(H,f,A,b,Aeq,beq,lb,ub)
   opts = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
   [z_op, fval] = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], opts);
   y = -Q_inv*(P*z_op + q);
   
   %% renormalize y:
   r = sqrt(y(1:3:end).^2+y(2:3:end).^2);
   r = repmat(r', [3, 1]);
   y = y./r(:);
   
   fv= 0.5*y'*A*y;
   hv= lambda*(norm(B(:,:,1)*y, 1)+norm(B(:,:,2)*y, 1)...
              +norm(B(:,:,3)*y, 1)+norm(B(:,:,4)*y, 1));
   
   fprintf('It %2d, fval=%f f=%f h=%f F=%f\n', it, fval, fv, hv, fv+hv);
   if (f_prev - (fv+hv) < 1e-5)
     % break;
   end
   f_prev = fv+hv;
end

