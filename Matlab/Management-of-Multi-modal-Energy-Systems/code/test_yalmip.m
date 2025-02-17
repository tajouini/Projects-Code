%cd c:\gurobi1002\win64\matlab
%addpath('c:\gurobi1002\win64\matlab')
% gurobi_setup
 if (~ isempty (getenv ('$SOFTWARE')))
    % we are on the cluster system
    addpath (getenv ('$SOFTWARE'));    % add the path for calling install_tools
    init_tools ();                  % initialize the tools
 else
        addpath('c:\gurobi1002\win64\matlab') % do whatever is required on your local file system
 end



sdpvar a x
Constraints = [a+1 <= x];
Objective = x^2;
ops = sdpsettings('solver','gurobi','verbose',0);
P = optimizer(Constraints,Objective,ops,a,x);
%Optimizer object with 1 inputs and 1 outputs. Solver: MOSEK-LP/QP