% simulation
 pol_index=[pol_indxH'; pol_indxL'];   
T_sim=5000;

a=[A_h; A_l];

%seed 
rng(1);
rand_nums=rand(T_sim,1);

%turn random numbers into values for A
%using transition matrix
A_sim=zeros(T_sim,1);
A_sim(1)=1;


for t=1:T_sim
   if  A_sim(t)==1 %the high state
             if  rand_nums(t)<0.977 
               A_sim(t+1) = 1; % keep on the high
             else 
                 A_sim(t+1) = 2; %go to the low state
             end
   elseif  rand_nums(t)<0.926 %the low state
               A_sim(t+1) = 2; % keep on the low
             else 
                 A_sim(t+1) =1;%go to the high state
             end
end   

  % A_sim(5001,:) = [];
%start with the arbitrary capital stok, then follow the policy function
%accordinly to simulated state in the current period
k_sim_index=zeros(T_sim,1);
k_sim_index(1)=5;

for t=1:T_sim
    k_sim_index(t+1)=pol_index(A_sim(t),k_sim_index(t));
end

A=zeros(T_sim,1);
k_sim=zeros(T_sim,1);


%filling up A with the actual value of A 
A=a(A_sim);
k_sim=k(k_sim_index);

%Then, I need to fill the information of y, given these results 

y_sim=(A'.*k_sim).^alpha;
%Get the standard deviation
sdy=std(y_sim)


