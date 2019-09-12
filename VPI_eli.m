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

cons_h = A_h*k_mat' .^ alpha + (1 - delta) * k_mat' - k_mat; 
cons_l = A_l*k_mat' .^ alpha + (1 - delta) * k_mat' - k_mat; 

ret_h = cons_h .^ (1 - sigma) / (1 - sigma); % return function high
ret_l = cons_l .^ (1 - sigma) / (1 - sigma); % return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret_h(cons_h < 0) = -Inf;
ret_l(cons_l < 0) = -Inf;

%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
%v_guess = zeros( 2, num_k);
v_guess = zeros(num_k, 2);
while dis > tol
    % compute the utility value for all possible combinations of k and k':
     [vfnL,pol_indxL]=max(ret_l + beta*repmat(v_guess*A(2,:)',1,num_k));
 [vfnH,pol_indxH]=max(ret_h + beta*repmat(v_guess*A(1,:)',1,num_k));

    %  value_matH = ret_h + beta * repmat(v_guess*A(1,:)', [1 num_k]);
    %   value_matL = ret_l + beta * repmat(v_guess*A(2,:)', [1 num_k]);
  % value_matH = ret_h + beta * repmat(A(1,:)*v_guess, [num_k 1]);
  % value_matL = ret_l + beta * repmat(A(2,:)*v_guess, [num_k 1]);
    % find the optimal k' for every k:
   % vfn = [vfnH'; vfnL'];
    vfn = [vfnH vfnL];
   
    vfn=[vfnH' vfnL'];
    % what is the distance between current guess and value function
    % based on the whole matrix
    dis = max(abs(vfn - v_guess));
       v_guess = vfn;
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
end

g = [k(pol_indxH)' k(pol_indxL)']; % policy function
plot(k,vfn,'Linewidth',1) 
xlabel('k') 
ylabel('V(k)')
title('Neoclasich Stochastic Value Functions Interaction')


plot(k,vfnH,k,vfnL)
figure
plot(k,g)


