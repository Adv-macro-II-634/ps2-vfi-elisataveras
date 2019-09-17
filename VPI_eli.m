clear 
close all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
A= [0.977 0.023; 0.074 0.926];
A_h=1.1;
A_l=.678;

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow

cons_h = A_h*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 
cons_l = A_l*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 

ret_h = cons_h .^ (1 - sigma) / (1 - sigma); % return function high
ret_l = cons_l .^ (1 - sigma) / (1 - sigma); % return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret_h(cons_h < 0) = -Inf;
ret_l(cons_l < 0) = -Inf;

%%%% Iteration
dis = 1; tol =1e-07; % tolerance for stopping 
v_guess = zeros(2, num_k);
while dis > tol
    % compute the utility value for all possible combinations of k and k':
     [vfnH,pol_indxH]=max(ret_h + beta*repmat(A(1,:)*v_guess,[num_k,1]),[],2);
    [vfnL,pol_indxL]=max(ret_l + beta*repmat(A(2,:)*v_guess,[num_k,1]),[],2);
    

   % vfn = [vfnH' vfnL'];
   vfn=[vfnH'; vfnL'];
    % what is the distance between current guess and value function
    % based on the whole matrix
    dis = max(abs(vfn - v_guess));
       v_guess = vfn;
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
end

g = [k(pol_indxH)' k(pol_indxL)']; % policy function
gH = k(pol_indxH); % policy function
gL= k(pol_indxL);
plot(k,vfnH','--',k,vfnL',':','Linewidth',1) 
xlabel('k') 
ylabel('V(k)')
title('Neoclasich Stochastic Value Functions Interaction')
legend({'A Hight','A Low'},'Location','southeast')
figure


%policy function
plot(k,k(pol_indxH)','--',k,k(pol_indxL)',':','Linewidth',1) 
xlabel('k') 
ylabel('V(k)')
title('VFI g(k) vs k')
legend({'A Hight','A Low'},'Location','southeast')



% PART 3 
kgrowthH=gH./k;
kgrowthL=gL./k;

plot(k,kgrowthH,'--',k,kgrowthL,':','Linewidth',1) 
xlabel('k') 
ylabel('V(k)')
title('VFI  results g(k)/k')
legend({'A Hight','A Low'},'Location','southeast')

%saving is f(k)-c

fH=A_h*k.^ alpha ;
fL=A_l*k.^ alpha ;

ch = A_h*k .^ alpha + (1 - delta) * k - gH; 
cl = A_l*k .^ alpha + (1 - delta) * k - gL; 

sH=(fH-ch)./k;
sL=(fL-cl)./k;


plot(k,sH,'--',k,sL,':','Linewidth',1) 
xlabel('k') 
ylabel('Saving')
title('VFI Saving per each Productivity Level')
legend({'A Hight','A Low'},'Location','southwest')


% PART 4
% First, I need to generate some random number for y, that has standard
% deviation of 0.018
% Then A=y/k.^alpha=
%now, in the data 


% PART 5

%using two loops over each 



clear 
close all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
A= [0.977 0.023; 0.074 0.926];
amat = [1.1 .678]'; % usng the matrix of amat

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit


polfun = zeros(num_k+1,3);
 v0 = zeros(num_k,1);
 dif = 10;
 its = 0;
maxits = 10;
 %%%% Iteration
dis = 1; tol =1e-07; % tolerance for stopping 
 while dis > tol & its < maxits
     for j = 1:3
     for i = 1:N
     k0 = kmat(i,1);
     a0 = amat(j,1);
     
     %%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow

c = a0*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 
ret = cons_h .^ (1 - sigma) / (1 - sigma); % return function high
ret_h(cons_h < 0) = -Inf;
ret_l(cons_l < 0) = -Inf;



     c = a0*k0.^alpha - k + (1-delta)*k0; % consumption
     
     k1 = fminbnd(valfunstoch,k_min,k_max);
     v1(i,j) = -vfstoch(k1);
     k11(i,j) = k1;
     end
     end
     %g = abs(v1?v0);
     dis = norm(v1-v0)
     v0 = v1;
     its = i
 end

figure
plot(kmat,v1,'?k','Linewidth',1)
xlabel('k')
ylabel('V(k)')

figure
plot(kmat,polfun,'?k','Linewidth',1)
xlabel('k')
ylabel('c')
