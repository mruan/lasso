clear;
dbstop if error;

%% Load Image:
img = imread('~/Dropbox/lasso/nov21/2lines.png');

%% Create the problem specific matrices
[L, C, c, D, N] = init_4way_sparse(img);

%% d'*x <=0 <=> c_i >= 0
%d  = zeros([3*N, 1]); d(3:3:end) = -1;
d  = zeros([N, 3*N]);
for i=1:N
   d(i, 3*i) = -1; 
end

%% Optimization specific parameters.
xi = 1e-05;  % noramally called lambda
mu = 1.0;

%% The problem is:
%% min (0.5*z'(A'/Q*A)z+(b'/Q*A+l)z  
%% where z = [u_{1-4} v w];
%% x_op = -Q\(Az+b)

%% these are independent of y
Q = L+mu*D; % Q = sparse(Q); 
Q_inv = 0.5*(eye(size(Q))/Q+ Q\eye(size(Q))); % To guarantee symmetry.

M = 14;
Aeq = []; beq = zeros(4*2*N,1);
for i=1:4*N
   Aeq = blkdiag(Aeq, [1 -1 0; 0 1 -1]); 
end
Aeq = [Aeq zeros(2*4*N, 2*N)];

%% l:
l = zeros([1, M*N]); l(end-N:end)=1;
lb = zeros([M*N,1]); ub = lb;

xi_range = 10.^(-2:0.2:0);
fv = zeros(size(xi_range));
hv = fv;
for k = 1:numel(xi_range)
    xi = xi_range(k);
    
%% upper and lower bounds
lb(1:(M-2)*N) = -xi; lb((M-1)*N+1:end) = -Inf; % -xi <= u and 0<=v
ub(1:(M-2)*N) =  xi; ub((M-2)*N+1:end) = +Inf; %   u <= xi 

%% The initial guess: y
y = zeros([3*N, 1]);
% y(1:3:end) = -1; y(3:3:end) =0.5;
y(1:3:end) = rand(N,1);
y(2:3:end) = sqrt(1-y(1:3:end).^2);
y(3:3:end) = rand(N, 1);

for it = 1:100
%% Make Y matrix
Y = zeros([3*N, N]);
for i=1:N
   Y(3*(i-1)+1: 3*(i-1)+2, i) = y(3*(i-1)+1: 3*(i-1)+2);
end

%% Only A and b depend on y
A = [C(:,:,1)', C(:,:,2)', C(:,:,3)', C(:,:,4)', d', Y];
A = sparse(A);
b = -mu*D*y;

%% let H=(A'/Q*A), f=(b'/Q*A+l)
H = A'*Q_inv*A; H = 0.5*(H+H'); % to guarantee symmetry
f = b'*Q_inv*A+l;

% x = quadprog(H,f,A,b,Aeq,beq,lb,ub)
opts = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
[z_op, fval] = quadprog(H, f, [], [], [], [], lb, ub, zeros(M*N, 1), opts);
y = -Q_inv*(A*z_op+b);

%% renormalize y:
r = sqrt(y(1:3:end).^2+y(2:3:end).^2);
r = repmat(r', [3, 1]);
y = y./r(:);

% fv= 0.5*y'*L*y;
% hv= xi*(norm(C(:,:,1)*y, 1)+norm(C(:,:,2)*y, 1)...
%                     +norm(C(:,:,3)*y, 1)+norm(C(:,:,4)*y, 1));
% fprintf('It %2d, fval=%f f=%f h=%f F=%f\n', it, fval, fv, hv, fv+hv);

end
fv(k) = 0.5*y'*L*y;
hv(k) = norm(C(:,:,1)*y, 1)+norm(C(:,:,2)*y, 1)+norm(C(:,:,3)*y, 1)+norm(C(:,:,4)*y, 1);
fprintf('Xi=%f, fv=%f, hv=%f\n', xi, fv(k), hv(k));
end
loglog(xi_range, fv, 'r', xi_range, hv, 'b');
% plot3(y(1:3:end), y(1:3:end), y(3:3:end), '.r');
% axis([-1 1 -1 1 0 2]); grid on; xlabel('a'); ylabel('b'); zlabel('z');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Aeq = []; beq = zeros(4*2*N,1);
% for i=1:4*N
%    Aeq = blkdiag(Aeq, [1 -1 0; 0 1 -1]); 
% end
% Aeq = [Aeq zeros(2*4*N, 2*N)];
