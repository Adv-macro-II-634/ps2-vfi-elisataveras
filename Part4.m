% simulation
T_sim=5000;

%seed 
rng(1);
rand_nums=rand(T_sim,1);

%turn random numbers into values for A
%using transition matrix
A_sim=zeros(T_sim,1);
A_sim(1)=1;

for t=i:T_sim
   if  A_sim(t)==1
             if  rand_nums(t)<0.977 < 0;
               A_sim(t+1) = 1;
             else 
                 A_sim(t+1) = 2;
             end
   else  
         if  rand_nums(t)<0.926 < 0;
               A_sim(t+1) = 1;
             else 
                 A_sim(t+1) = 2;
             end
   end   



%start with the arbitrary capital stok, then follow the policy function
%accordinly to simulated state in the current period
k_sim_index=zeros(T_sim,1);
k_sim_index(1)=5;

for t=i:T_sim
    k_sim(t+1)=pol_index(A_sim(t),k_sim_index(t));
end

%How to simulate state A 



