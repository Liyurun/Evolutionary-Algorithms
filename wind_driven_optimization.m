function position = wind_driven_optimization(obj, problem)

opt = struct('max_iter', 10,...
                'pop_size' , 20,... 
                'npar',  5,...		 
                'maxit',  500,...		 
                'RT',  3,...			 
                'g',  0.2,...			 
                'alp',  0.4,...		 
                'c',  0.4,...			 
                'maxV',  0.3,...			 
                'disp',1,...
                'problem',problem);
        wdo = struct('name','WDO');
        wdo = opt;
        wdo = start(wdo);
        
            if wdo.disp
                disp(['WDO  ' 'started ...']);
                disp('Initializing population.');
            end       
        
        wdo = initialize(wdo);
        
        % Iterations (Main Loop)
                for it = 1:wdo.max_iter
                    
                    % Set Iteration Counter
                    wdo.iter = it;
                    
                    % Iteration Started
                    wdo = iterate(wdo);
                    
                    % Update Histories
                    wdo.position_history(it).position = wdo.best_sol.position;
                    wdo.nfe_history(it) = wdo.nfe;
                    wdo.best_obj_value_history(it) = wdo.best_sol.obj_value;
                    
          
                    % Display Iteration Information
                    if wdo.disp
                        disp(['Iteration ' num2str(wdo.iter) ...
                              ': Best de. Value = ' num2str(wdo.best_sol.obj_value) ...
                              ', position = ' num2str(wdo.best_sol.position)]);
                    end
                    
                    position =wdo.best_sol.position;

                end
                
end








        % Initialization
        function this = initialize(this)
            
            % Create Initial Population (Not Sorted)
            sorted = false;
                        this.best_sol.obj_value = inf;
            this = init_pop(this,sorted);
            

            
        end
        
        % Iterations
        function this = iterate(this)
            
          
            % Decision Vector Size
            var_size = [1 this.problem.dim];

            % Create New Population
            newpop = repmat(this.empty_individual, this.pop_size, 1);
            [globalpres,indx] = min([this.pop(:).obj_value]);
            globalpos = this.pop(indx).position;
            % Update the velocity:
            for i=1:this.pop_size
            % choose random dimensions:
            a = randperm(this.problem.dim);        			
            % choose velocity based on random dimension:
                velot(i,:) = this.vel(i,a);				
                this.vel(i,:) = (1-this.alp)*this.vel(i,:)-(this.g*this.pop(i).position)+ ...
                        abs(1-1/i)*((globalpos-this.pop(i).position).*this.RT)+ ...
                        (this.c*velot(i,:)/i);
            end
            
            
            for j = 1:size(this.pop,1)
                this.pop(j).position = this.pop(j).position + this.vel(j);
                this.pop(j).obj_value = this.problem.func(this.pop(j).position);
            end
            
            
        % Finding best particle in population
    	[minpres,indx] = min([this.pop.obj_value]);
    	
        
    	%----------------------------------------------------
    	% Rank the air parcels:
    	[sorted_pres rank_ind] = sort([this.pop.obj_value]);
    	% Sort the air parcels position, velocity and pressure:
        this.vel = this.vel(rank_ind,:);
        for j = 1:size(this.pop,1)
    	this.pop(j).position = this.pop(rank_ind).position;
    	
    	this.pop(j).obj_value = sorted_pres(j);  
        end
        
        
    	% Updating the global best:
    	better = minpres < this.best_sol.obj_value;
    	if better
            this.best_sol.obj_value = minpres;
            this.best_sol.position = this.pop(indx).position;           	% min location for this iteration
        end
end
            
 
 
    
    
    
    
 

        
        
        
        
        
        
        
        
        
        
        
        
        % Reset the Algorithm
function this = start(this)
            
this.empty_individual.position = [];
this.empty_individual.obj_value = [];
this.empty_individual.solution = [];
this.pop= [];
this.best_sol= [];
this.iter= 0;
this.nfe= 0;
this.run_time= 0;
this.avg_eval_time= 0;
this.nfe_history= [];
this.best_obj_value_history= [];
this.must_stop= 0;
% reset

            this.iter = 0;
            this.nfe = 0;
            this.nfe_history = nan(1, this.max_iter);
            this.best_obj_value_history = nan(1, this.max_iter);
            this.run_time = 0;
            this.avg_eval_time = 0;
            this.must_stop = false;
end
        
        % Decode and Evaluate a Coded Solution
        function [z, sol] = decode_and_eval(this, xhat)
            
            % Increment NFE
            this.nfe = this.nfe + 1;
            
            % Decode and Evaluate

                [z, sol] = this.problem.func(xhat);

    
        end
        
        % Evaluate a Single Solution or Population
        function pop = eval(this, pop)
            for i = 1:numel(pop)
                pop(i).position = ypea_clip(pop(i).position, 0, 1);
                [pop(i).obj_value, pop(i).solution] = this.decode_and_eval(pop(i).position);
            end
        end
        
        % Create a New Individual
        function ind = new_individual(this, x)
            
            if ~exist('x', 'var') || isempty(x)
                x = rand([1 this.problem.dim]);
            end
            
            ind = this.empty_individual;
            ind.position = ypea_clip(x, 0, 1);
            [ind.obj_value, ind.solution] =  decode_and_eval(this,x);
        end
        
        % Initialize Population
        function this = init_pop(this, sorted)
            
            % Check for Sorted flag (Default is false, not sorted)
            if ~exist('sorted', 'var') || isempty(sorted)
                sorted = false;
            end
            
            % Initialize Population Array
            this.pop = repmat(this.empty_individual, this.pop_size, 1);
            this.vel = this.maxV * 2 * (rand(this.pop_size,this.problem.dim)-0.5);  

            % Initialize Best Solution to the Worst Possible Solution
            this.best_sol = this.empty_individual;
            this.best_sol.obj_value = this.problem.worst_value;
            
            % Generate New Individuals
            for i = 1:this.pop_size
                
                % Generate New Solution
                this.pop(i) = new_individual(this);
                
                % Compare to the Best Solution Ever Found
                if ~sorted && is_better(this,this.pop(i), this.best_sol)
                    this.best_sol = this.pop(i);
                end
                
            end
            
            % Sort the Population if it is needed
            if sorted
                this.pop = sort_population(this,this.pop);
                this.best_sol = this.pop(1);
            end
            
        end
        
        % Sort Population
        function [pop, sort_order, obj_values] = sort_population(this, pop)
            
 
            % Sort the Objective Values Vector
            [obj_values, sort_order] = sort([pop.obj_value] );
            
            % Sort (Re-order) Population
            pop = pop(sort_order);
            
        end
        
        % Sort and Select the Population
        function pop = sort_and_select(this, pop)
            
            % Sort Population
            pop = sort_population(this,pop);
            
            % Set the Population Size Limit
            n = min(this.pop_size, numel(pop));
            pop = pop(1:n);
            
        end
        
        % Check if a solution is better than other
        function b = is_better(this, x1, x2)
            b = x1.obj_value<x2.obj_value;
        end
        
        % Get Best Memebr of Population
        function pop_best = get_population_best(this, pop)
            pop_best = pop(1);
            for i = 2:numel(pop)
                if this.is_better(pop(i), pop_best)
                    pop_best = pop(i);
                end
            end
        end
        
        % Get Positions Matrix of Population
        function pos = get_positions(this, pop)
            pos = reshape([pop.position], this.problem.dim, [])';
        end
        
        % Get Objective Values of Population
        function v = get_objective_values(~, pop)
            v = [pop.obj_value];
        end
        
        % Calculate Selection Probabilities
        function p = get_selection_probs(this, values, selection_pressure)
            
            % Check if Selection Pressure specified
            if ~exist('selection_pressure', 'var')
                selection_pressure = 1;
            end
            
            % Change selection pressure sign, according to problem type

                alpha = -selection_pressure;
            
            % Calculate Selection Probabilities
            p = exp(alpha*values);
            p = ypea_normalize_probs(p);
            
        end






