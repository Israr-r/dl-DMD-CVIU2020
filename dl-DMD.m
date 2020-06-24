% example MATLAB script for loading neurofinder data
%
% for more info see:
%
% - http://neurofinder.codeneuro.org
% - https://github.com/codeneuro/neurofinder
%
% requires one package from the matlab file exchange
%
% - jsonlab
% - http://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files-in-matlab-octave
%/Users/israr/Downloads/Fluo-N2DH-SIM+/01/t000.tif
%/Users/israr/Downloads/Fluo-N2DH-SIM+/01_GT/SEG/man_seg000.tif
% load the images
%files = dir('images/*.tiff');
clear all
clc
%%
imgs = []
%/Users/israr/Downloads/Fluo-N2DH-SIM+/01_GT/TRA


cd '/Users/israr/Desktop/Dl-DMD-Coefficients/dlDMDvsStdDMD';
load vid_72_90
cd '/Users/israr/Desktop/VISAPP2020/Visapp_2020_code';
%vid = vid(:,:,3501:3700);
zt = vid;
%zt = vid+0.8*rand(size(zt));


Xi = zt;
%Xi = Xi(:,:,1:2:end);

nof1 = size(Xi,3); 
nof = size(Xi,3); % # of frames
p_size = 8; % patch size
X = [];

for i = 1:nof1 % for DMD not for dl-DMD
     temp = Xi(:,:,i);
    X(:,i) = temp(:)';     %  making column matrix of input video
    %count = count + 1 ;  
end

m  = size(Xi,1);
n  = size(Xi,2);


t = linspace(0,16*pi,size(Xi,3));  % time points
dt = t(2) - t(1);  


%% Dictionary Learning 

img = [Xi(:,:,1) Xi(:,:,60) Xi(:,:,120)] % select random patches from these frames

%p_size = 8;
col1 = im2col(img,[p_size p_size],'sliding');

%ind = round(1 + (size(col1,2)-1).*rand(10000,1));

No_patches_train = 5000;
col = col1(:,randperm(No_patches_train));
col = col/max(col(:));

Natoms = 64;
Beta_d1 = zeros(Natoms,size(col,2));
D1 = col1(:,randperm(Natoms));

%D = zeros(size(D));
err= [];
fn = [];

D = rand(p_size*p_size,Natoms);


for i = 1:100


Beta_d1 = beta_calc1_dict_learn(col,D,Beta_d1);
D = my_dict_learn(col,D,Beta_d1);
% for j = 1:Natoms
%    D(:,j) = D(:,j)/max(D(:,j));
%     
% end

%D = D/max(D(:));

fn = norm((col-D*Beta_d1));
err = [err fn];
i

if(fn<0.0001)
   break; 
end
end
figure(121),plot(err)
Dict = D;

%% Estimating coefficients
ytemp = [];
temp = [];
y = [];
vecmean = [];
for i = 1 : nof % Aligning patches vertically   
   % ytemp(:,:,i) = im2col(Xi(:,:,i),[p_size p_size],'sliding');
    ytemp(:,:,i) = im2col(Xi(:,:,i),[p_size p_size],'sliding');
  
    temp = ytemp(:,:,i);
    y(:,i) = temp(:); % size = (64 x #  of patches) x (# of frames (N))   
    %vecmean(:,i) = repmat(mean(temp(:)),[d 1]); 
end

clear ytemp

%% Coefficient matrix estimation
No_patches_m = size(im2col(Xi(:,:,1),[p_size p_size],'sliding'),2);

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




%% sliding reconstruction check

%X1p = zeros(size(X0));% Approximated Signals
%2p = zeros(size(X0));
%% For overlap patches 
% Aligned verticle video frames
B1 = Beta1;%(:,1:end-1);
B2 = Beta2;%(:,2:end);

%X1p = zeros(size(X0));% Approximated Signals
%2p = zeros(size(X0));
%% For overlap patches 
% Aligned verticle video frames
for i = 1: No_patches_m
   X1p((i-1)*(p_size.^2)+1:i*(p_size.^2),:) = D*Beta1((i-1)*(Natoms)+1:i*(Natoms),:);
    
end
 
 for i = 1: No_patches_m
    X2p((i-1)*(p_size.^2)+1:i*(p_size.^2),:) = D*Beta2((i-1)*(Natoms)+1:i*(Natoms),:);    
 end

% for distinct columns

% chk=[];
% X3p = [];
% for i = 1: No_patches_m
% X3p((i-1)*(p_size.^2)+1:i*(p_size.^2),:) = D*Phi((i-1)*(Natoms)+1:i*(Natoms),:);
% 
% end


%% DMD 

%X1p = Beta1; %To apply DMD on coefficients
%X2p = Beta2;
% noise_beta = 0*rand(size(Beta1));
% X1p = Beta1 + noise_beta ;
% X2p = Beta2 + noise_beta ;
X1p = Beta1;
X2p = Beta2;

[U, S, V] = svd(X1p, 'econ');
r = nof-10; % number of frames -1
U_r = U(:, 1:r); % truncate to rank-r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);

% [idx] = find(diag(S)>1e-8);
% r = nof-1;
% U_r = U(:, idx); % truncate to rank-r
% S_r = S(idx, idx);
% V_r = V(:, idx);

Atilde = U_r' * X2p * V_r / S_r; % low-rank dynamics
%Atilde = U_r' * B2 * (pinv(S_r') * V_r')';

[W_r, D] = eig(Atilde);
Phi = X2p * V_r / S_r * W_r; % DMD modes
%Phi = B2 * (pinv(S_r') * V_r')' * W_r; % DMD modes

lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues
abs(omega)

%%  Background subtraction
bg = find(abs(omega)<1e-2); % Threshold for background separattion
%bg = find(max(abs(omega(:))));
%bg = find(abs(omega)>20);
%bg = find(abs(omega)<1e-2);
%fg = setdiff(idx, bg);
%fg = find(abs(omega)<20);
fg = setdiff(1:r, bg);

omega_fg = omega(fg); % foreground
Phi_fg = Phi(:,fg); % DMD foreground modes

omega_bg = omega(bg); % background
Phi_bg = Phi(:,bg); % DMD background mode
%%  Background subtraction
bg = find(abs(omega)<1e-2); % Threshold for background separattion
%bg = find(abs(omega)<1e-2);
%fg = setdiff(idx, bg);
fg = setdiff(1:size(omega,1), bg);

omega_fg = omega(fg); % foreground
Phi_fg = Phi(:,fg); % DMD foreground modes

omega_bg = omega(bg); % background
Phi_bg = Phi(:,bg); % DMD background mode




%% Full reconstruction

D = Dict;

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
  %%A = A+ ones(size(A,1),1)*vecmean(i,:);
  tempbg = my_im2coldmd(A,Xi(:,:,i),2,1); % sigma  
  %tempbg = col2im(A,[p_size p_size],[m n], 'distinct'); % sigma  
 % Ibg(:,:,i) = real(tempbg/max(tempbg(:)));
  I_full(:,:,i) = real(tempbg/max(tempbg(:)));

end
implay(real(I_full))


%% Background

b = Phi_bg \ X1p(:, 1); % amplitude

X_bg = zeros(numel(omega_bg), length(t));
for tt = 1:length(t),
    X_bg(:, tt) = b .* exp(omega_bg .* t(tt));
end;
X_bg = real(Phi_bg * X_bg);
X_bg = X_bg(1:end, :);
%% 
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
  tempbg = col2im(A,[p_size p_size],[m n], 'distinct'); % for 
 % Ibg(:,:,i) = real(tempbg/max(tempbg(:)));
  Ibg(:,:,i) = tempbg;%real(tempbg/max(tempbg(:)));
i
end
implay(real(Ibg))


%% Foreground 

b = Phi_fg \X1p(:, 1); % amplitudes
%b = pinv(Phi_fg)*B1(:, 1);
X_fg = zeros(numel(omega_fg), length(t));
for tt = 1:length(t),
    X_fg(:, tt) = b .* exp(omega_fg .* t(tt));
end;
X_fg = real(Phi_fg * X_fg);
X_fg = X_fg(1:end, :);
%save('Ibg009.mat','Ibg','-v7.3');
%% Foreground reconstruction 

Xbgreconstruct = zeros(size(X0,1),nof);
 for i = 1: No_patches_m
    Xfgreconstruct((i-1)*(p_size.^2)+1:i*(p_size.^2),:) = D*X_fg((i-1)*(Natoms)+1:i*(Natoms),:);
    i 
 end

 % for distinct patches
 for i = 1:nof
  A = reshape((Xfgreconstruct(:,i)),p_size^2,No_patches_m);
  %%A = A+ ones(size(A,1),1)*vecmean(i,:);
  %tempbg = my_im2coldmd(A,Xi(:,:,i),p_size,1); % sigma  
  tempbg = col2im(A,[p_size p_size],[m n], 'distinct'); % sigma  
 % Ibg(:,:,i) = real(tempbg/max(tempbg(:)));
  Ifg(:,:,i) = real(tempbg/max(tempbg(:)));

end
implay(real(-Ifg))


 
 
 
 % for sliding patches

for i = 1:nof
  A = reshape((Xfgreconstruct(:,i)),p_size^2,No_patches_m);
  %%A = A+ ones(size(A,1),1)*vecmean(i,:);
  tempbg = my_im2coldmd(A,Xi(:,:,i),p_size,1); % sigma  
  %tempbg = col2im(A,[p_size p_size],[m n], 'distinct'); % sigma  
 % Ibg(:,:,i) = real(tempbg/max(tempbg(:)));
  Ifg(:,:,i) = real(tempbg/max(tempbg(:)));

end
implay(real(-Ifg))


%% making video 
v = VideoWriter('Background_dmd_Ibg.avi');  % video name
v.FrameRate = 15; % frame rate
open(v)
for i = 1:nof
    writeVideo(v,real(im_norm(real(Ibg(:,:,i))))); % Ibg = background, Ifg = Foreground
end
close(v)


%% Full reconstruction

b = Phi_bg \ X1p(:, 1); % amplitude

X_bg = zeros(numel(omega_bg), length(t));
for tt = 1:length(t),
    X_bg(:, tt) = b .* exp(omega_bg .* t(tt));
end;
X_bg = real(Phi_bg * X_bg);
X_bg = X_bg(1:end, :);


Phi_n = [];
for i = 1 : length(t)
   Phi_n(:,:,i) = reshape(X_bg(:,i),[400,400]);
    
end

figure,imshow(my_map(permute(abs(Phi_n(:,:,1:end)),[3,1,2])),[])
colormap hot
implay(abs(permute(Phi_n(:,end:-1:1,:),[2,1,3]))/max(real(Phi_n(:))))
colormap hot





bfull = Phi \ X1p(:,1);
X_full = zeros(numel(omega), length(t));
for tt = 1:length(t),
    X_full(:, tt) = bfull .* exp(omega .* t(tt));
end;
time_dynamics_full = real(Phi * X_full);

Phi_n = [];
for i = 1 : r
   Phi_n(:,:,i) = reshape(time_dynamics_full(:,i),[400,400]);
    
end

figure,imshow(my_map(permute(abs(Phi_n(:,:,1:2)),[3,1,2])),[])
colormap hot
implay(abs(Phi_n)/max(real(Phi_n(:))))
colormap hot






%% SVD truncation

[U, S, V] = svd(X1p, 'econ');
r = 2; % number of frames -1
U_r = U(:, 1:r); % truncate to rank-r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
rec = U_r*S_r*V_r';

svd_n = [];

for i = 1 : size(rec,2)
   svd_n(:,:,i) = reshape(rec(:,i),[400,400]);    
end

figure,imshow(my_map(permute(abs(svd_n(:,:,1:255)),[3,1,2])),[])
colormap hot
%% ------ Compressive sensing ------------
n = 400;
K = 5; % degree of sparsity
T = 2; % duration of integration
dt = 0.01; % time-step
r = 15; % number of PCA modes to keep
saveFLAG = 0
%%
M = 300;  % number of measurements

% projType = 1; % uniform random projection
projType = 2; % Gaussian random projection
% projType = 3; % Single pixel measurement

C = zeros(M,n*n);
Theta = zeros(M,n*n);
xmeas= zeros(n,n);
for i=1:M
    xmeas = 0*xmeas;
    if(projType==1)
        xmeas = rand(n,n);
    elseif(projType==2)
        xmeas = randn(n,n);
    elseif(projType==3)
        xmeas(ceil(n*rand),ceil(n*rand)) = 1;
    end
    C(i,:) = reshape(xmeas,n*n,1);
    Theta(i,:) = reshape((ifft2(xmeas)),1,n*n);
end

%% Project data
XDAT = X;
XD1 = X(:,1:end-1);
XD2 = X(:,2:end);
YDAT = C*XDAT;
if(saveFLAG)
    save([filename,'_DATAY.mat']);
end

YD1 = YDAT(:,1:end-1);
YD2 = YDAT(:,2:end);
% step 1
[U,S,V] = svd(YD1,0);
% step 2
r=10; % for two mode and K=5
Sinv = S(1:r,1:r)^(-1);
Atilde = U(:,1:r)'*YD2*V(:,1:r)*Sinv(1:r,1:r);
% step 3
[W,D] = eig(Atilde);
% step 4
PhiY = YD2*V(:,1:r)*Sinv*W;
PhiXr = XD2*V(:,1:r)*Sinv*W;  % reconstructed
%% Reconstruct Modes with Compressive Sensing.
for i=1:K
    i
    y = PhiY(:,2*i-1);    
%     cvx_begin;
%         variable ahat(n*n,1);
%         minimize( norm (ahat,1) );
%         subject to
%         Theta*ahat == y;
%     cvx_end;
    ahat = cosamp(Theta,y,10,10^-5,100);
    Phi(:,2*i-1) = reshape(real(ifft2(reshape(ahat,n,n))),n*n,1);
    Phi(:,2*i) = reshape(imag(ifft2(reshape(ahat,n,n))),n*n,1);
end

Phi_c = [];
for i = 1 : size(Phi,2)
   Phi_c(:,:,i) = reshape(Phi(:,i),[400,400]);
    
end
implay(abs(Phi_c/max(Phi_c(:))))
figure,imshow(my_map(permute(abs(Phi_c(:,:,1:10)),[3,1,2])),[])
colormap hot


%% plot C
figure, colormap bone
Cimg = sum(C,1);
Ximg = XDAT(:,1).*Cimg';
imagesc(reshape(Ximg,n,n));
figure, colormap bone;
surf(X,Y,Z,reshape(Cimg,n,n));
view(10,32)
axis equal, axis off
drawnow

%% ------



lambda = diag(D); % discrete-time eigenvalues
omega = log(lambda)/dt; % continuous-time eigenvalues
abs(omega)

%%  Background subtraction
%bg = find(abs(omega)<1e-1); % Threshold for background separattion
%bg = find(max(abs(omega(:))));
bg = find(abs(omega)>20);
%bg = find(abs(omega)<1e-2);
%fg = setdiff(idx, bg);
fg = find(abs(omega)<20);
%fg = setdiff(1:r, bg);

omega_fg = omega(fg); % foreground
Phi_fg = Phi(:,fg); % DMD foreground modes

omega_bg = omega(bg); % background
Phi_bg = Phi(:,bg); % DMD background mode


%% Full reconstruction
bfull = Phi \ X1p(:,1);
X_full = zeros(numel(omega), length(t));
for tt = 1:length(t),
    X_full(:, tt) = bfull .* exp(omega .* t(tt));
end;
time_dynamics_full = real(Phi * X_full);


% optimize Phi
%Phi_tilde = dict_calc1(X,X_full);
%beta_calc1(X,D);

%X_full = X_full(1:end, :);
%full_vid = Phi_tilde*X_full;

%% Play full video
full_vid = real(reshape(time_dynamics_full,[m,n,nof]));
implay(full_vid/full_vid(:))

[ssim_v,psnr_v, mse_v] = calc_graph_comparison(full_vid,vid)
%% Compute DMD Background Solution
b = Phi_bg \ X1p(:, 1); % amplitude

X_bg = zeros(numel(omega_bg), length(t));
for tt = 1:length(t),
    X_bg(:, tt) = b .* exp(omega_bg .* t(tt));
end;
X_bg = real(Phi_bg * X_bg);
X_bg = X_bg(1:end, :);
%% 

for i = 1:nof
  A = reshape(X_bg(:,i),[m,n]);
  %A = A+ ones(size(A,1),1)*vecmean(i,:);
  %tempbg = my_im2coldmd(A,Xi(:,:,i),8,1); % sigma  
 % Ibg(:,:,i) = real(tempbg/max(tempbg(:)));
  Ibg(:,:,i) = A;%real(tempbg/max(tempbg(:)));

end
implay(Ibg/max(Ibg(:)))

%% Compute DMD Foreground Solution
b = Phi_fg \X1p(:, 1); % amplitudes
%b = pinv(Phi_fg)*B1(:, 1);
X_fg = zeros(numel(omega_fg), length(t));
for tt = 1:length(t),
    X_fg(:, tt) = b .* exp(omega_fg .* t(tt));
end;
X_fg = real(Phi_fg * X_fg);
X_fg = X_fg(1:end, :);
%save('Ibg009.mat','Ibg','-v7.3');
%% Foreground reconstruction 

for i = 1:nof
  A = reshape(X_fg(:,i),[m,n]); %
 % A = A+ ones(size(A,1),1)*vecmean(i,:);
  % temp= my_im2coldmd(A,Xi(:,:,i),8,1); % sigma  
   Ifg(:,:,i) = (A)/max(A(:));%real(temp/max(temp(:)));
end
implay(Ifg/max(Ifg(:))) 

  At = reshape((y(:,1)),p_size^2,No_patches_m);
  %%A = A+ ones(size(A,1),1)*vecmean(i,:);
  %tempbg = my_im2coldmd(A,Xi(:,:,i),8,1); % sigma  
  tempbg = col2im(A,[p_size p_size],[m n], 'distinct'); % sigma  


%% -------------------------------------
% load the regions (training data only)
regions = loadjson('regions/regions.json');
masks = zeros(x, y);

for i = 1:length(regions)
	if isstruct(regions)
		coords = regions(i).coordinates;
	elseif iscell(regions)
		coords = regions{i}.coordinates;
	end
	masks(sub2ind([x, y], coords(:,1), coords(:,2))) = 1;
end

% show the outputs
figure();
subplot(1,2,1);
imagesc(mean(imgs,3)); colormap('gray'); axis image off;
subplot(1,2,2);
imagesc(masks); colormap('gray'); axis image off;
