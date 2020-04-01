function run()
% =============================================================================
% Optimized variables
%       theta = {x_1, x_2, ...}
%
% There are five types of algorithms that are integrated into this classdef, ranging
% from deterministic, heuristic algorithms to Bayesian Inference:
%       - Particle Swarm Optimizatio (PSO)
%       - Differential Evolution (DE)
%       - Markov chain Monte Carlo (MCMC) (Jacobian matrix might be needed)
%       - Metropolis Adjusted Differential Evolution (MADE)
%       - Metropolis Adjusted Langevin Algorithm (MALA)
%           defined on the Riemann geometry, and combined with parallel tempering
%           Jacobian matrix must be needed
% =============================================================================


    % There are four optimization algorithms availabe in this programme
    optimization_method = struct('Particle_Swarm_Optimization',0,...
                                 'Differential_Evolution',0,...
                                 'Genetic_Algorithm',0,...
                                 'Covariance_Matrix_Adaptation_Evolution_Strategy',0,...
                                 'Simulated_Annealing',1,...
                                 'Metropolis_Adjusted_Differential_Evolution',0,...
                                 'Parallel_Riemann_Metropolis_Adjusted_Langevin',0,...
                                 'Markov_Chain_Monte_Carlo',0,...
                                 'Deterministic_algorithm_fmincon',0);


    % The initial boundary of parameters: In the format of [x^1_min x^1_max; ...]
    opt.paramBound = [-32.768 -12.768; 
                      32.768 12.768];
    problem = struct('func',@(x) func(x),...
                     'upbound',opt.paramBound(1,:),...
                     'lobound',opt.paramBound(2,:),...
                     'dim',size(opt.paramBound(1,:),2),...
                     'worst_value',Inf);
      
    if optimization_method.Particle_Swarm_Optimization

        OptAlgorithms.Particle_Swarm_Optimization(opt, problem);

    elseif optimization_method.Differential_Evolution

        OptAlgorithms.Differential_Evolution(opt, problem);

    elseif  optimization_method.Genetic_Algorithm

        OptAlgorithms.Genetic_Algorithm(opt, problem);
        
    elseif  optimization_method.Covariance_Matrix_Adaptation_Evolution_Strategy
        
        OptAlgorithms.Covariance_Matrix_Adaptation_Evolution_Strategy(opt, problem);
       
    elseif  optimization_method.Simulated_Annealing
        
        OptAlgorithms.Simulated_Annealing(opt, problem);
        %%
    elseif isfield(optimization_method, 'Markov_Chain_Monte_Carlo')

        OptAlgorithms.Markov_Chain_Monte_Carlo(opt, problem);

    elseif isfield(optimization_method, 'Metropolis_Adjusted_Differential_Evolution') 

        OptAlgorithms.Metropolis_Adjusted_Differential_Evolution(opt, problem);



    elseif isfield(optimization_method, 'Deterministic_algorithm_fmincon') ...
            && optimization_method.Deterministic_algorithm_fmincon

        % The guessing point for the deterministic method 
        initParams = [1 1];

        % Check the consistence of the initial boundary condition and the parameter amount
        OptAlgorithms.checkOptDimension(opt, length(initParams));

        loBound = opt.paramBound(:,1);
        upBound = opt.paramBound(:,2);

        options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter',...
            'TolX',1e-6,'TolCon',1e-6,'TolFun',1e-6,'MaxIter',500);

        try
            [xValue, yValue, exitflag, output, ~, grad] = fmincon( @objectiveFunc, ...
                initParams, [],[],[],[], loBound, upBound, [], options);
        catch exception
            disp('Errors in the MATLAB build-in optimizer: fmincon. \n Please check your input parameters and run again. \n');
            disp('The message from fmincon: %s \n', exception.message);
        end

        fprintf('----------------  Minimum: %10.3g  ---------------- \n', yValue);
        fprintf('%10.3g | ', xValue);
        fprintf('\n------------------------------------------------------ \n');

    else

        error('The method you selected is not provided in this programme \n');

    end

end
% =============================================================================
%              The MATLAB library for optimization case studies
% 
%      Copyright © 2015-2016: Qiaole He
% 
%      Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
% 
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
