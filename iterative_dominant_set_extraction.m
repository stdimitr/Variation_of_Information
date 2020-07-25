function [groups,no_groups,cost_function,f_ini,Q]=iterative_dominant_set_extraction(A)
%
%  [groups,no_groups,cost_function_values,f_ini]=iterative_dominant_set_extraction(A)
%  
%  Subtractive Clustering based on the function listed below
%
%  groups contains class labels : e.g. [1 2 1 1 3 1 0 1 1] 
%   ---> zero denotes spurious nodes in the graph
%   different graph-components have been extracted sequentially 
%
%   cost_function --> tabulates the corresponding 'cluster-quality' for each component
%
%    f_ini ---> cluster quality for the whole graph taken as single component
%
%   It is based on the function
%  [sel_list,rest_list,ordered_list,memberships]=dominant_set_extraction(A)
%                  A is a similarity (weighted adjacency) matri 
%
%  The algorithm of which has been adopted from
%  PAMI, vol.29(1),2007,pp.167 ---> Dominant Sets & Pairwise Clustering
%

%http://users.auth.gr/laskaris/
%http://users.auth.gr/stdimitr/


clear {A,groups,rem_list,f_ini,cost_function,zn,z,f_current,no_groups,sel_list,rest_list,group_new} 



A=A-diag(diag(A)); % to enfroce no self loop for each node
[N,N]=size(A);


groups(1,N)=0;
rem_list(1:N)=0;
no_groups=0;
cost_function=0;
f_ini=0;
f_current=0;
zn(1:N)=0;
z(1:N)=0;

groups=zeros(1,N);
rem_list=[1:N];
i=1;
cost_function=[];
zn=ones(1,N);
zn=zn/sum(zn);
f_ini=zn*A*zn'; 
f_current=f_ini+0.5*f_ini;

k=0;

while (f_current-f_ini)*(1-isempty(rem_list))>0        %f_current >= f_ini
   k=k+1;
[sel_list,rest_list]=dominant_set_extraction(A(rem_list,rem_list),k);
group_new=rem_list(sel_list); 
groups(group_new)=i; 
z=zeros(1,N);
z(group_new)=1;
z=z/sum(z); 
f_current=z*A*z';
cost_function(i)=f_current;
rem_list=setdiff(rem_list,group_new);
i=i+1;
end
no_groups=i-1;


K=sum(A);                               %degree
m=sum(K);                               %number of edges (each undirected edge is counted twice)
B=A-(K.'*K)/m;                    %modularity matrix
                           %community indices

s=groups;                      %compute modularity
Q=~(s-s.').*B/m;
Q=sum(Q(:));



