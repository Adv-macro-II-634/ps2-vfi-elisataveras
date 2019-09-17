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
amat = [A_h A_l]'; 

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
v_guessL = zeros(1, num_k);
v_guessH = zeros(1, num_k);
%v_guess = zeros(num_k, 2);
while dis > tol
    %generating the continuation for each possible state 
    HH= beta*repmat(0.977*v_guessH,[num_k,1]);
    HL=beta*repmat(0.023*v_guessH,[num_k,1]);
    LH= beta*repmat(0.074*v_guessL,[num_k,1]);
    LL= beta*repmat(0.926*v_guessL,[num_k,1]);
    [vfnH,pol_indxH]=max(ret_h+HH+HL,[],2);
    [vfnL,pol_indxL]=max(ret_l+LH+LL,[],2); 
   
    vfnHigh=vfnH' ;
    vfnLow=vfnL' ;
    vfn = [vfnH' vfnL'];
   %vfn=[vfnH'; vfnL'];
   % vfn=[vfnH'];
    % what is the distance between current guess and value function
    % based on the whole matrix
    dH=abs(vfnHigh - v_guessH); 
    dL=abs(vfnLow - v_guessL); 
    dis =  max(dH,dL)
       v_guessH = vfnHigh;
       v_guessL = vfnLow;
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
end

gH=k(pol_indxH);
gL=k(pol_indxL);

plot(k,vfnHigh,'--',k,vfnLow,':','Linewidth',1) 
xlabel('k') 
ylabel('V(k)')
title('Neoclasich Stochastic Value Functions Interaction')

figure

plot(k,gH,k,gL,'Linewidth',1) 
xlabel('k') 
ylabel('V(k)')
title('VFI g(k) vs k')


% PART 3 
kgrowthH=gH./k;
kgrowthL=gL./k;

plot(k,kgrowthH,k,kgrowthL,'Linewidth',1) 
xlabel('k') 
ylabel('V(k)')
title('VFI g(k)/k VFI')


%saving is f(k)-c

fH=A_h*k.^ alpha ;
fL=A_l*k.^ alpha ;

ch = A_h*k .^ alpha + (1 - delta) * k - gH; 
cl = A_l*k .^ alpha + (1 - delta) * k - gH; 

sH2=fH-ch;
sL2=fL-cl;



plot(k,sH2,k,sL2,'Linewidth',1) 
xlabel('k') 
ylabel('V(k)')
title('VFI g(k)-k/k VFI')



% or is t k'-k

sH=(gH-k)./k;
sL=(gL-k)./k;


% or saving rate z=k'/f(k)

zH=gH./fH;
zL=gL./fL;


sH=zH./k;
sL=zL./k;


plot(k,sH,k,sL,'Linewidth',1) 
xlabel('k') 
ylabel('V(k)')
title('VFI g(k)-k/k VFI')


% PART 4 SIMULATION OF A
% only 1 random sequence for A
% the long run probability are 0.237113402 for Low and 0.762886598
% for High
% it is still true that A_h=(1-0.237113402A_l)/.762886598
%now, I need to simulate a random seuqnece for A that I use to get A_h
%This needs to match the output std of 1.8% 
%V(A)=V(Y/k^alpha)
% creat a random number from a normal distribution, whtib ;
% sign the variance of two variabel
% v(y)=V(A)V(k^alpha)+Var(A)(E(k^alpha))^2+Var(k^alpha)(E(A))^2
% V(A)=v(y)-Var(k^alpha))/(V(k^alpha)+(E(k^alpha))^2)

%solve for stadand deviation of A 
A_g = zeros(1, num_k);
eqn=std(A_g.*k.^alpha)==0.018;
sA=solve(eqn,A_g)
vy= (0.018)^2; 
vk= var(k.^alpha);
ek= mean(k.^alpha)^2;
vA=(vy-vk)/(vk+ek)
sA=sqrt(vA);

% PART 5



%%%% Iteration


