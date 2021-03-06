function position =biogeography_based_optimization(obj, problem)

    opt = struct('max_iter', 10,...
                'pop_size' , 20,...
                'keep_rate' , 0.5,...
                'alpha' , 0.9,...
                'mutation_prob' , 0.1,...
                'mutation_step_size' , 0.2,...
                'mutation_step_size_damp' , 0.98,...
                'disp',1,...
                'problem',problem);
        bbo = struct('name','BBO');
        bbo = opt;
        bbo = start(bbo);
        
            if bbo.disp
                disp(['BBO  ' 'started ...']);
                disp('Initializing population.');
            end       
        
        bbo = initialize(bbo);
        
        % Iterations (Main Loop)
                for it = 1:bbo.max_iter
                    
                    % Set Iteration Counter
                    bbo.iter = it;
                    
                    % Iteration Started
                    bbo = iterate(bbo);
                    
                    % Update Histories
                    bbo.position_history(it).position = bbo.best_sol.position;
                    bbo.nfe_history(it) = bbo.nfe;
                    bbo.best_obj_value_history(it) = bbo.best_sol.obj_value;
                    
          
                    % Display Iteration Information
                    if bbo.disp
                        disp(['Iteration ' num2str(bbo.iter) ...
                              ': Best de. Value = ' num2str(bbo.best_sol.obj_value) ...
                              ', position = ' num2str(bbo.best_sol.position)]);
                    end
                    
                    position =bbo.best_sol.position;

                end
end


        % Initialization
        function this = initialize(this)
            
            % Create Initial Population (Sorted)
            sorted = true;
            this = init_pop(this,sorted);
            
            % Set Keep Count
            this.params.keep_count = round(this.keep_rate * this.pop_size);
            
            % Set Newly Created Solutions Count
            this.params.new_count = this.pop_size - this.params.keep_count;
            
            % Emmigration Rates
            this.params.mu = linspace(1, 0, this.pop_size);
            
            % Immigration Rates
            this.params.lambda = 1 - this.params.mu;
            
            % Initial Value of Mutation Step Size
            this.params.sigma = this.mutation_step_size;
            
        end
        
        % Iterations
        function this = iterate(this)
            
            % Decision Vector Size
            var_count = this.problem.dim;
            
            % Create New Population
            newpop = this.pop;
            for i = 1:this.pop_size
                
                % Generate New Solution
                xnew = newpop(i).position;
                for k = 1:var_count
                    
                    % Migration
                    if rand <= this.params.lambda(i)
                        
                        % Emmigration Probabilities
                        EP = this.params.mu;
                        EP(i) = 0;
                        EP = EP/sum(EP);
                        
                        % Select Source Habitat
                        j = roulette_wheel_selection(EP);
                        
                        % Migration
                        xnew(k) = this.pop(i).position(k) ...
                            + this.alpha*(this.pop(j).position(k) - this.pop(i).position(k));
                        
                    end

                    % Mutation
                    if rand <= this.mutation_prob
                        xnew(k) = xnew(k) + this.params.sigma*randn;
                    end
                    
                end
                
                % Create New Solution
                newpop(i) =  new_individual(this,xnew);
                
            end
            
            % Sort, Select and Merge
            keep_count = this.params.keep_count;
            new_count = this.params.new_count;
            newpop =  sort_population(this,newpop);
            this.pop =  sort_and_select(this,[this.pop(1:keep_count); newpop(1:new_count)]);
            
            % Update Best Solution Ever Found
            this.best_sol = this.pop(1);
            
            % Damp Mutation Step Size
            this.params.sigma = this.mutation_step_size_damp * this.params.sigma;
                        
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
                pop(i).position = clip(pop(i).position, 0, 1);
                [pop(i).obj_value, pop(i).solution] = this.decode_and_eval(pop(i).position);
            end
        end
        
        % Create a New Individual
        function ind = new_individual(this, x)
            
            if ~exist('x', 'var') || isempty(x)
                x = rand([1 this.problem.dim]);
            end
            
            ind = this.empty_individual;
            ind.position = clip(x, 0, 1);
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
            p = normalize_probs(p);
            
        end
        
        
        function y = clip(x, lb, ub)
    % Clips the inputs, and ensures the lower and upper bounds.
    
    if ~exist('lb', 'var')
        lb = 0;
    end
    if ~exist('ub', 'var')
        ub = 1;
    end
    
    y = min(max(x, lb), ub);
    
end


function p = normalize_probs(p)
    % Normalize Probabilities
    
    p(p<0) = 0;
    
    p = p/sum(p);
    
    p(isnan(p)) = 0;
    
    if all(p==0)
        p(:) = 1/numel(p);
    end
    
end



function L = roulette_wheel_selection(P, count, replacement)
    % Performs Roulette Wheel Selection    
    
    if ~exist('count', 'var')
        count = 1;
    end

    if ~exist('replacement','var')
        replacement = false;
    end    
    
    if ~replacement
        count = min(count, numel(P));
    end
    
    C = cumsum(P);
    S = sum(P);
    
    L = zeros(count, 1);
    for i = 1:count
        L(i) = find(rand()*S <= C, 1, 'first');
        if ~replacement
            P(L(i)) = 0;
            C = cumsum(P);
            S = sum(P);
        end
    end
    

end