
classdef OptAlgorithms < handle
% =============================================================================
% This is the class of the functions of optimization algorithms.
%
% =============================================================================   

    properties (Constant)
        FUNC        = @func;   % the objective function
        swarm       = 2;               % the population amount
        sample      = 10;             % the loop count for evolution
        dataPoint   = 20;              % the amount of data observation
        DE_strategy = 5;                % the strategy of DE kernel
    end

%   DE
    methods (Static = true, Access = 'public')

        function position = Differential_Evolution(obj, problem)
            position = differential_evolution(obj, problem);

    end % DE
    end

%   PSO    
    methods (Static = true, Access = 'public')

        function position = Particle_Swarm_Optimization(obj, problem)

            position = particle_swarm_optimization(obj, problem);
            
        end
    end


%   GA    
    methods (Static = true, Access = 'public')

        function position = Genetic_Algorithm(obj, problem)

            position = genetic_algorithm(obj, problem);
            
        end
    end
    
    %   CMAES   
    methods (Static = true, Access = 'public')

        function position = Covariance_Matrix_Adaptation_Evolution_Strategy(obj, problem)

            position =covariance_matrix_adaptation_evolution_strategy(obj, problem);
            
        end
    end
    % SA 
    methods (Static = true, Access = 'public')

        function position = Simulated_Annealing(obj, problem)

            position =simulated_annealing(obj, problem);
            
        end
    end
     % ABC
    methods (Static = true, Access = 'public')

        function position = Artificial_Bee_Colony(obj, problem)

            position =artificial_bee_colony(obj, problem);
            
        end
    end
         % ACO
    methods (Static = true, Access = 'public')

        function position =    Ant_Colony_Optimization(obj, problem)

            position = ant_colony_optimization(obj, problem);
            
        end
    end
             % BA
    methods (Static = true, Access = 'public')

        function position =    Bees_Algorithm(obj, problem)

            position = bees_algorithm(obj, problem);
            
        end
    end
      % BBO
    methods (Static = true, Access = 'public')

        function position = Biogeography_Based_Optimization(obj, problem)

            position = biogeography_based_optimization(obj, problem);
            
        end
    end
      % FA
    methods (Static = true, Access = 'public')

        function position = Firefly_Algorithm(obj, problem)

            position = firefly_algorithm(obj, problem);
            
        end
    end

          % HS
    methods (Static = true, Access = 'public')

        function position = Harmony_Search(obj, problem)

            position = harmony_search(obj, problem);
            
        end
    end
    
              % ICA
    methods (Static = true, Access = 'public')

        function position = Imperialist_Competitive_Algorithm(obj, problem)

            position = imperialist_competitive_algorithm(obj, problem);
            
        end
    end
    
    
   % IWO
    methods (Static = true, Access = 'public')

        function position = Invasive_Weed_Optimization(obj, problem)

            position = invasive_weed_optimization(obj, problem);
            
        end
    end
    
    
    % TLBO
    methods (Static = true, Access = 'public')

        function position = Teaching_Learning_Based_Optimization(obj, problem)

            position = teaching_learning_based_optimization(obj, problem);
            
        end
    end
    
    
  %WOA
    methods (Static = true, Access = 'public')

        function position = Whale_Optimization_Algorithm(obj, problem)

            position = whale_optimization_algorithm(obj, problem);
            
        end
    end
    
    
%WDO
    methods (Static = true, Access = 'public')

        function position = Wind_Driven_Optimization(obj, problem)

            position = wind_driven_optimization(obj, problem);
            
        end
    end
    
    
%FSA
    methods (Static = true, Access = 'public')

        function position = Future_Search_Algorithm(obj, problem)

            position = future_search_algorithm(obj, problem);
            
        end
    end
    
end
