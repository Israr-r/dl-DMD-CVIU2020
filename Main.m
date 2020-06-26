
clear all
clc
%%

load vid

Xi = vid;

nof1 = size(Xi,3); 
nof = size(Xi,3); % # of frames
p_size = 8; % patch size
X = [];

for i = 1:nof1 % for DMD not for dl-DMD
     temp = Xi(:,:,i);
     X(:,i) = temp(:)';  % making column matrix of input video
   
end

[m,n,nof] = size(Xi)

t = linspace(0,16*pi,size(Xi,3));  % time points
dt = t(2) - t(1);  


%% Dictionary Learning 
a = 1;
b = nof;
r = round((b-a).*rand(3,1) + a); % Generate random 
img = [Xi(:,:,r(1)) Xi(:,:,r(2)) Xi(:,:,r(3))]; % select random patches from these frames

col1 = im2col(img,[p_size p_size],'sliding'); % Generate patches from image 

No_patches_train = 5000; % Number of patches to train dictionary
col = col1(:,randperm(No_patches_train)); % Create random patches
col = col/max(col(:));

Natoms = 64; % No. of dictionary atoms
Beta_d1 = zeros(Natoms,size(col,2));
%D1 = col1(:,randperm(Natoms)); 

err= [];
fn = [];

D = rand(p_size*p_size,Natoms); % Start with randomly generated dictionary

% Dictionary learning
for i = 1:50 % Number of iterations to learn dictionary

Beta_d1 = beta_calc1_dict_learn(col,D,Beta_d1);
D = my_dict_learn(col,D,Beta_d1);
fn = norm((col-D*Beta_d1));
err = [err fn];

if(fn<0.0001)
   break; 
end
end
figure(2),plot(err), title('Loss for dictionary learning')
Dict = D; % Dictionary 

%% Estimating coefficients
ytemp = [];
temp = [];
y = [];
vecmean = [];
for i = 1 : nof % Aligning patches vertically   
   % ytemp(:,:,i) = im2col(Xi(:,:,i),[p_size p_size],'sliding'); 
     ytemp(:,:,i) = im2col(Xi(:,:,i),[p_size p_size],'distinct');  
     temp = ytemp(:,:,i);
     y(:,i) = temp(:); % size = (64 x #  of patches) x (# of frames (N))   
end

clear ytemp

%% Coefficient matrix estimation
No_patches_m = size(im2col(Xi(:,:,1),[p_size p_size],'distinct'), 2);

X0 = y(:,1:end-1); 
Beta1=zeros(Natoms*No_patches_m,nof-1);

for i = 1:No_patches_m
      %Beta1((i-1)*(Natoms)+1:i*Natoms,:) =OMPerr(D,X0((i-1)*(p_size.^2)+1:i*(p_size.^2),:),10e-20);
     Beta1((i-1)*(Natoms)+1:i*Natoms,:) = beta_calc1(X0((i-1)*(p_size.^2)+1:i*(p_size.^2),:),D);
     i
end
%clear X0

X1 = y(:,2:end); % Frames -> 2: number of frames (nof)
Beta2=zeros(Natoms*No_patches_m,nof-1);

for i = 1:No_patches_m
    %Beta2((i-1)*(Natoms)+1:i*Natoms,:) =OMPerr(D,X1((i-1)*(p_size.^2)+1:i*(p_size.^2),:),10e-20);
    Beta2((i-1)*(Natoms)+1:i*Natoms,:) = beta_calc1(X1((i-1)*(p_size.^2)+1:i*(p_size.^2),:),D);
    i
end
%clear X1



%% DMD

X1p = Beta1;
X2p = Beta2;

[U, S, V] = svd(X1p, 'econ');
r = nof-5; % number of frames -1
U_r = U(:, 1:r); % truncate to rank-r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);

Atilde = U_r' * X2p * V_r / S_r; % low-rank dynamics

[W_r, D] = eig(Atilde);
Phi = X2p * V_r / S_r * W_r; % DMD modes

lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues

%%  Background subtraction
bg = find(abs(omega)<1e-2); % Threshold for background separation
fg = setdiff(1:r, bg);
omega_fg = omega(fg); % foreground
Phi_fg = Phi(:,fg); % DMD foreground modes
omega_bg = omega(bg); % background
Phi_bg = Phi(:,bg); % DMD background mode



%% Full Reconstruction of Video Frames

D = Dict; % dictionary 

bfull = Phi \ X1p(:,1);
X_full = zeros(numel(omega), length(t));
for tt = 1:length(t),
    X_full(:, tt) = bfull .* exp(omega .* t(tt));
end;
time_dynamics = real(Phi*X_full);

X0 = y(:,1:end-1); % Frames -> 1: (number of frames -1)

 X_full = zeros(size(X0,1),nof);
 for i = 1: No_patches_m
    X_full((i-1)*(p_size.^2)+1:i*(p_size.^2),:) = D*time_dynamics((i-1)*(Natoms)+1:i*(Natoms),:);
    i 
 end

for i = 1:nof
  A = reshape((X_full(:,i)),p_size^2,No_patches_m);
  %tempbg = my_im2coldmd(A,Xi(:,:,i),2,1); % sigma  
  tempbg = col2im(A,[p_size p_size],[m n], 'distinct'); % sigma  
  %Ibg(:,:,i) = real(tempbg/max(tempbg(:)));
  I_full(:,:,i) = real(tempbg/max(tempbg(:)));
end
implay(real(I_full))


%% Background Reconstruction

b = Phi_bg \ X1p(:, 1); % amplitude
X_bg = zeros(numel(omega_bg), length(t));
for tt = 1:length(t),
    X_bg(:, tt) = b .* exp(omega_bg .* t(tt));
end;
X_bg = real(Phi_bg * X_bg);
X_bg = X_bg(1:end, :);

X0 = y(:,1:end-1); % Frames -> 1: (number of frames -1)


 Xbgreconstruct = zeros(size(X0,1),nof);
 for i = 1: No_patches_m
    Xbgreconstruct((i-1)*(p_size.^2)+1:i*(p_size.^2),:) = D*X_bg((i-1)*(Natoms)+1:i*(Natoms),:);
    i 
 end


for i = 1:nof
  A = reshape((Xbgreconstruct(:,i)),p_size^2,No_patches_m);
  %%A = A+ ones(size(A,1),1)*vecmean(i,:);
  %tempbg = my_im2coldmd(A,Xi(:,:,i),8,1); % sigma  
  tempbg = col2im(A,[p_size p_size],[m n], 'distinct'); 
 % Ibg(:,:,i) = real(tempbg/max(tempbg(:)));
  Ibg(:,:,i) = tempbg;%real(tempbg/max(tempbg(:)));
i
end
implay(real(Ibg)) % Reconstructed Background


%% Foreground reconstruction 

b = Phi_fg \X1p(:, 1); % amplitudes
%b = pinv(Phi_fg)*B1(:, 1);
X_fg = zeros(numel(omega_fg), length(t));
for tt = 1:length(t),
    X_fg(:, tt) = b .* exp(omega_fg .* t(tt));
end;
X_fg = real(Phi_fg * X_fg);
X_fg = X_fg(1:end, :);

Xbgreconstruct = zeros(size(X0,1),nof);
 for i = 1: No_patches_m
    Xfgreconstruct((i-1)*(p_size.^2)+1:i*(p_size.^2),:) = D*X_fg((i-1)*(Natoms)+1:i*(Natoms),:);
    i 
 end

 for i = 1:nof
  A = reshape((Xfgreconstruct(:,i)),p_size^2,No_patches_m);
  %tempbg = my_im2coldmd(A,Xi(:,:,i),p_size,1); % sigma  
  tempbg = col2im(A,[p_size p_size],[m n], 'distinct'); % sigma  
  %Ibg(:,:,i) = real(tempbg/max(tempbg(:)));
  Ifg(:,:,i) = real(tempbg/max(tempbg(:)));
end
implay(-Ifg) % Reconstructed Foreground
