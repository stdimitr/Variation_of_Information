

%A small example of how to run vi matlab function
%If you employ vi function please refer to :
%(1)Dimitriadis SI, Laskaris NA, Del Rio-Portilla Y, Koudounis GC. 
% Characterizing dynamic functional connectivity across sleep stages from EEG. 
% Brain Topography Volume 22, Number 2 / September, 2009 p.119-133. 
% or
% (2)Dimitriadis SI, Laskaris NA, Tsirka V, Vourkas V, Micheloyannis S.
% An EEG study of brain connectivity dynamics at the resting state. 
% Nonlinear Dynamics, Physchology and Life Sciences.

%1)create 2 random connectivity graphs
nodes=30;
graph1=rand(30,30);
graph2=rand(30,30);

for i=1:nodes
    graph1(i,i)=0; 
    graph2(i,i)=0;
end



for k=1:nodes
   for l=(k+1):nodes
       graph1(l,k)=graph1(k,l);
       graph2(l,k)=graph2(k,l);
   end
end

%2)run dominant sets algorithm for each graph seperately
[groups1,no_groups,cost_function,f_ini,Q]=iterative_dominant_set_extraction(graph1);

[groups2,no_groups,cost_function,f_ini,Q]=iterative_dominant_set_extraction(graph2);

%3)run vi to get the distance between two clusterings
[VI_value,NVI, adjvi] = vi(groups1',groups2');
