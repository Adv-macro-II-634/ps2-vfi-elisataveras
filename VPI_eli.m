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
%v_guess = zeros(num_k, 2);
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

% PART 4
% First, I need to generate some random number for y, that has standard
% deviation of 0.018
% Then A=y/k.^alpha=
%now, in the data 


% PART 5


k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

cons_h = A_h*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 
cons_l = A_l*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 


dis = 1; tol = 1e-06; % tolerance for stopping 
v_guess = zeros(1, num_k);
while dis > tol
    for j = 1:2 %space of possible states
    for i = 1:num_k % space of k values 
      K0 = k_mat(i,1);
      A0 = amat(j,1);
     end
    end
end





