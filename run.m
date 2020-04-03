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
                                 'Simulated_Annealing',0,...
                                 'Artificial_Bee_Colony',0,...
                                 'Ant_Colony_Optimization',0,...
                                 'Bees_Algorithm',0,...
                                 'Biogeography_Based_Optimization',0,...
                                 'Firefly_Algorithm',0,...
                                 'Harmony_Search',0,...
                                 'Imperialist_Competitive_Algorithm',0,...
                                 'Invasive_Weed_Optimization',0,...
                                 'Teaching_Learning_Based_Optimization',0,...
                                 'Whale_Optimization_Algorithm',0,...
                                 'Wind_Driven_Optimization',0,...
                                 'Future_Search_Algorithm',1);


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
        
    elseif  optimization_method.Artificial_Bee_Colony
        
        OptAlgorithms.Artificial_Bee_Colony(opt, problem);
    
    elseif  optimization_method.Ant_Colony_Optimization
        
        OptAlgorithms.Ant_Colony_Optimization(opt, problem);
        
    elseif  optimization_method.Bees_Algorithm
        
        OptAlgorithms.Bees_Algorithm(opt, problem);
        
    elseif  optimization_method.Biogeography_Based_Optimization
        
        OptAlgorithms.Biogeography_Based_Optimization(opt, problem);
    elseif  optimization_method.Firefly_Algorithm
        
        OptAlgorithms.Firefly_Algorithm(opt, problem);
    elseif  optimization_method.Harmony_Search
        
        OptAlgorithms.Harmony_Search(opt, problem);                
      
    elseif  optimization_method.Imperialist_Competitive_Algorithm
        
        OptAlgorithms.Imperialist_Competitive_Algorithm(opt, problem);      
        
    elseif  optimization_method.Invasive_Weed_Optimization
        
        OptAlgorithms.Invasive_Weed_Optimization(opt, problem);  
        
    elseif  optimization_method.Teaching_Learning_Based_Optimization
        
        OptAlgorithms.Imperialist_Competitive_Algorithm(opt, problem);  
        
    elseif  optimization_method.Whale_Optimization_Algorithm
        
        OptAlgorithms.Whale_Optimization_Algorithm(opt, problem);  
        
        
    elseif  optimization_method.Wind_Driven_Optimization
        
        OptAlgorithms.Wind_Driven_Optimization(opt, problem);  
        
    elseif  optimization_method.Future_Search_Algorithm
        
        OptAlgorithms.Future_Search_Algorithm(opt, problem); 
        
        
        
    elseif  optimization_method.Markov_Chain_Monte_Carlo
        
        OptAlgorithms.Markov_Chain_Monte_Carlo(opt, problem);  
        
        
        
        
    elseif isfield(optimization_method, 'Markov_Chain_Monte_Carlo')

        OptAlgorithms.Markov_Chain_Monte_Carlo(opt, problem);

    elseif isfield(optimization_method, 'Metropolis_Adjusted_Differential_Evolution') 

        OptAlgorithms.Metropolis_Adjusted_Differential_Evolution(opt, problem);
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
