function x1 = beta_calc1_dict_learn(y,H,x)


lambda = 1e-5;
alpha = 0;
alphak = 1e5;


T = lambda/(2*alphak);

tk = 1;
tk1 = 1;
   err =[]; 
   
fn = inf;
fp = inf;

yk1 = x;
x1 = x;
x2 = [];
x0 = x;
for k = 1:500%    number of iterations
    
gk = -2*H'*(y-H*yk1);
x1 = yk1 - gk/alphak;
x1 = soft(x1, T);
tk1 = (1 + sqrt(1+4*tk.^2))/2;
yk1 = x1 + ((tk-1)/(tk1)).*(x1-x0);
tk = tk1;

fn = norm((y-H*x1)) + lambda*norm(x1,1);
err = [err fn];

 e = norm(x1 - x0)/numel(x1);
 if e < 10e-500
     %break;
  end


x0 = x1;
fp = fn;

if mod(k,500)==0
    plot(err), title('Beta Coefficeint')
    pause(0.5)
end

end

