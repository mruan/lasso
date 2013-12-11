function [A, B, b, D, N, beta_idx_lin] = init_4way_sparse(img)
[H, W, Depth] = size(img);
if(Depth==3) % if rgb image -> grayscale image
  img = rgb2gray(img);
end

%% Convert to BW binary image
bw_img = img < 128;

beta_idx_lin = find(bw_img);  % we are interested in these pixels...
N = numel(beta_idx_lin);      % dimension of the problem
[beta_idx_u, beta_idx_v] = ind2sub([H, W], beta_idx_lin);

A = zeros(3*N);
D = zeros(3*N);
B = zeros([3*N, 3*N, 4]);
b = zeros([  N, 3*N, 4]);

% helper function sub2idx
sub2idx =@(y, x) find(beta_idx_lin == sub2ind([H, W], y, x));

% Find the normalizing scale
sc = 1/max([W, H]); 
for i=1:N
  u = beta_idx_u(i);
  v = beta_idx_v(i);
  
  k = (i-1)*3+1:i*3;
  A(k, k) = [u*sc v*sc 1]'*[u*sc v*sc 1];
  D(k, k) = diag([1 1 0]);
  
  % Find all neighbors of pixel indexed at i:
  if(u>1 && bw_img(u-1,v))
    t = sub2idx(u-1, v); r = (t-1)*3+1:t*3;
    B(k, k, 1) = diag([1,1,1]);
    B(k, r, 1) = diag([-1,-1,-1]);
    b(i, k, 1) =  1;
    b(i, r, 1) = -1;
  end
  
  if(v>1 && bw_img(u,v-1))
    t = sub2idx(u, v-1); r = (t-1)*3+1:t*3;
    B(k, k, 2) = diag([1,1,1]);
    B(k, r, 2) = diag([-1,-1,-1]);
    b(i, k, 2) =  1;
    b(i, r, 2) = -1;
  end
  
  if(u<H && bw_img(u+1,v))
    t = sub2idx(u+1, v); r = (t-1)*3+1:t*3;
    B(k, k, 3) = diag([1,1,1]);
    B(k, r, 3) = diag([-1,-1,-1]);
    b(i, k, 3) =  1;
    b(i, r, 3) = -1;
  end
  
  if(v<W && bw_img(u,v+1))
    t = sub2idx(u, v+1); r = (t-1)*3+1:t*3;
    B(k, k, 4) = diag([1,1,1]);
    B(k, r, 4) = diag([-1,-1,-1]);
    b(i, k, 4) =  1;
    b(i, r, 4) = -1;
  end
end

end

