function x1 = my_dict_learn(y,x,H)


%%

%x = zeros(size(y*H'));

%lambda = 0.0000012;
lambda = 1e-5;
alpha = 0;
alphak = 1e5;
T = lambda/(2*alphak);
tk = 1;
tk1 = 1;
err =[];    
% load H
% load x
% load y

fn = inf;
fp = inf;

yk1 = x; 

%yk1 = zeros(size(x));
%yk1 = x;
x1 = x;%zeros(size(yk1));
x2 = x;%zeros(size(x1));
x0 = x;%zeros(size(x1));
for k = 1:100%    
 % phi(k)
 gk = -2*(y-yk1*H)*H'; % for dictionary 
 x1 = yk1 - gk/alphak;
 x1 = soft(x1, T);
 tk1 = (1 + sqrt(1+4*tk.^2))/2;
 yk1 = x1 + ((tk-1)/(tk1)).*(x1-x0);
 tk = tk1;
% x1 = yk1 - gk/alphak;
% x1 = soft(x1, T);
% tk1 = (1 + sqrt(1+4*tk.^2))/2;
% yk1 = x1 + ((tk-1)/(tk1)).*(x1-x0);
% tk = tk1;

%fn = norm((y-x1*H),'fro') + lambda*norm(H,1);
fn = norm((y-x1*H),'fro') + lambda*norm(H,1);
err = [err fn];

 e = norm(x1 - x0)/numel(x1);
 if e < 10e-15
    % break;
  end

% if fn>fp
%    break; 
% end
x0 = x1;
fp = fn;
%x1 = x2;

%alpha = ((x1 - x0)'*(gk - gkm1))/((x1 - x0)'*(x1 - x0)); 
%gkm1 = gk;

%x = wthresh(x,'s',T); 

% gk = - H' * (y - H*x);r
% x = soft(x - gk/alphak, T);
% x2=x;
k
if mod(k,1)==0
plot(err),title('Dictionary')
pause(0.5)
check = 1;
end
end

%figure,plot(err)
