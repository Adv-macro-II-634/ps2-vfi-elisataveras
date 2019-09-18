% simulation
 pol_index=[p_indexH; p_indexL];
        
T_sim=5000;

%seed 
rng(1);
rand_nums=rand(T_sim,1);

%turn random numbers into values for A
%using transition matrix
A_sim=zeros(T_sim,1);
A_sim(1)=1;

for t=1:T_sim
   if  A_sim(t)==1
             if  rand_nums(t)<0.977 
               A_sim(t+1) = 1;
             else 
                 A_sim(t+1) = 2;
             end
   elseif  rand_nums(t)<0.926
               A_sim(t+1) = 2;
             else 
                 A_sim(t+1) = 1;
             end
   end   

%start with the arbitrary capital stok, then follow the policy function
%accordinly to simulated state in the current period
k_sim_index=zeros(T_sim,1);
k_sim_index(1)=5;

for t=1:T_sim
    k_sim_index(t+1)=pol_index(A_sim(t),k_sim_index(t));
end

A=zeros(T_sim,1);
k_sim=zeros(T_sim,1);

%Now, fill the data with the actual values of the simulation A_h=1.1;
%A_l=.678;
if A_sim==1
    A=1.1;
else 
  A=.678;
end
  
  k_sim=K(k_sim_index);
k_sim(:,5001) = [];
%Then, I need to fill the information of y, given these results 

ksimT=k_sim';
At=A';

y_sim=(At.*ksimT).^alpha;
%Get the standard deviation
sdy=std(y_sim)

