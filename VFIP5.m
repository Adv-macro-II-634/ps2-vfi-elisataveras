close all

alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
A_h=1.1;
A_l=.678;
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k
N = 1000;
tol = 1e-05; % tolerance for convergence of the value functions
maxits=100;
k_bar = delta ^ (- 1 / (1 - alpha));
k_max = k_bar;

K = linspace(0, k_max, N+1);
K(1) = []; % an N x 1 column vector for the space of capital

V_guessH = zeros(1,N); % a row vector for the initial guess of the value function
V_guessL = zeros(1,N); % a row vector for the initial guess of the value function
% V_guessH =  0.9 * K' ; % alternative guess

dis = 1; % measure of distance of value functions to each other
n = 0; % count number of iteration steps
 its = 0;
tic
while dis > tol & its < maxits
    
    V_newH = zeros(N,1); % this will be the updated function (result of applying operator to 'V_guessH')
    policy_fnH = zeros(N,1);
    V_newL = zeros(N,1); % this will be the updated function (result of applying operator to 'V_guessH')
    policy_fnL = zeros(N,1);
 
    
    for ik = 1:N % loop over values for k
        
        k_today = K(ik);
        
        valueH_ik = zeros( N, 1 );
        valueL_ik = zeros( N, 1 );
        for ii = 1:N % loop over k tomorrow
            
       k_prime = K(ii);
        cH=A_h*k_today^ alpha + (1 - delta) * k_today - k_prime;
        cL=A_l*k_today^ alpha + (1 - delta) * k_today - k_prime;
       valueH_ik(ii) = (cH).^ (1 - sigma) / (1 - sigma) + (beta*0.977*V_guessH(ii))+(beta*0.023*V_guessL(ii));
       valueL_ik(ii) = (cL).^ (1 - sigma) / (1 - sigma) + (beta*0.074*V_guessH(ii))+(beta*0.926*V_guessL(ii));
             
            % need to make sure that resource constraint isn't violated! If
            % that happens, we can just set the value for such a kprime to
            % something very negative to make sure it will never be chosen
            if cH < 0
                valueH_ik(ii) = -realmax;
            end
            if cL < 0
                valueL_ik(ii) = -realmax;
            end
            
        end
        
        [best_valueH, indexH] = max(valueH_ik);
        [best_valueL, indexL] = max(valueL_ik);
        
        V_newH(ik) = best_valueH;
        policy_fnH(ik) = K(indexH);
        
          V_newL(ik) = best_valueL;
        policy_fnL(ik) = K(indexL);
        
    end
    
    % compute distance b/w V~n and V~n+1
    dH = max( abs( V_guessH - V_newH' ) );
    dL = max( abs( V_guessL - V_newL' ) );
    dis = max(dH,dL);
    
    % update guess
    V_guessH = V_newH';
     V_guessL = V_newL';
       its = i  ;
end
toc