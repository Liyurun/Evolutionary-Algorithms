function position = artificial_bee_colony(obj, problem)

    opt = struct('max_iter', 10,...
                'pop_size' , 20,...
                'onlooker_count',22,...
                'max_acceleration',0.5,...
                'disp',1,...
                'problem',problem);
        abc = struct('name','ABC');
        abc = opt;
        abc = start(abc);
        
            if abc.disp
                disp(['ABC  ' 'started ...']);
                disp('Initializing population.');
            end       
        
        abc = initialize(abc);
        
        % Iterations (Main Loop)
                for it = 1:abc.max_iter
                    
                    % Set Iteration Counter
                    abc.iter = it;
                    
                    % Iteration Started
                    abc = iterate(abc);
                    
                    % Update Histories
                    abc.position_history(it).position = abc.best_sol.position;
                    abc.nfe_history(it) = abc.nfe;
                    abc.best_obj_value_history(it) = abc.best_sol.obj_value;
                    
          
                    % Display Iteration Information
                    if abc.disp
                        disp(['Iteration ' num2str(abc.iter) ...
                              ': Best de. Value = ' num2str(abc.best_sol.obj_value) ...
                              ', position = ' num2str(abc.best_sol.position)]);
                    end
                    
                    position =abc.best_sol.position;

                end
                

        
        
        
end


function abc = start(abc)
 
abc.empty_individual.position = [];
abc.empty_individual.obj_value = [];
abc.empty_individual.solution = [];
abc.pop= [];
abc.best_sol= [];
abc.iter= 0;
abc.nfe= 0;
abc.run_time= 0;
abc.avg_eval_time= 0;
abc.nfe_history= [];
abc.best_obj_value_history= [];
abc.must_stop= 0;
% reset

            abc.iter = 0;
            abc.nfe = 0;
            abc.nfe_history = nan(1, abc.max_iter);
            abc.best_obj_value_history = nan(1, abc.max_iter);
            abc.run_time = 0;
            abc.avg_eval_time = 0;
            abc.must_stop = false;

end

        % Initialization
        function abc = initialize(abc)
            
            % Create Initial Population (Not Sorted)
            sorted = false;
            abc = init_pop(abc,sorted);
            
            % Abandonment Limit Parameter (Trial Limit)
            abc.params.L = round(0.6 * abc.problem.dim * abc.pop_size);
            
            % Abandonment Counter
            abc.params.C = zeros(abc.pop_size, 1);
            
        end
        
        % Iterations
        function abc = iterate(abc)
            
            % Decision Vector Size
            var_size = [1 abc.problem.dim];
            
            % Recruited Bees
            for i = 1:abc.pop_size
                
                % Choose k randomly, not equal to i
                K = [1:i-1 i+1:abc.pop_size];
                k = K(randi([1 numel(K)]));
                
                % Define Acceleration Coeff.
                phi = abc.max_acceleration*uniform_rand(-1 , 1, var_size);

                % New Bee Position
                xnew = abc.pop(i).position ...
                     + phi.*(abc.pop(i).position - abc.pop(k).position);
                newsol = new_individual(abc,xnew);
                
                % Comparision
                if is_better(abc,newsol, abc.pop(i))
                    abc.pop(i) = newsol;
                else
                    abc.params.C(i) = abc.params.C(i)+1;
                end

            end

            % Calculate Fitness Values and Selection Probabilities
            obj_values = get_objective_values(abc,abc.pop);
            P = get_selection_probs(abc,obj_values);
            
            % Onlooker Bees
            for m = 1:abc.onlooker_count

                % Select Source Site
                i = roulette_wheel_selection(P);
                
                % Choose k randomly, not equal to i
                K = [1:i-1 i+1:abc.pop_size];
                k = K(randi([1 numel(K)]));
                
                % Define Acceleration Coeff.
                phi = abc.max_acceleration*uniform_rand(-1 , 1, var_size);

                % New Bee Position
                xnew = abc.pop(i).position ...
                     + phi.*(abc.pop(i).position - abc.pop(k).position);
                newsol = new_individual(abc,xnew);
                
                % Comparision
                if is_better(abc,newsol, abc.pop(i))
                    abc.pop(i) = newsol;
                else
                    abc.params.C(i) = abc.params.C(i)+1;
                end
                
            end

            % Scout Bees
            for i = 1:abc.pop_size
                if abc.params.C(i) >= abc.params.L
                    abc.pop(i) = new_individual(abc);
                    abc.params.C(i) = 0;
                end
            end
            
            % Update Best Solution Ever Found
            for i = 1:abc.pop_size
                if is_better(abc,abc.pop(i), abc.best_sol)
                    abc.best_sol = abc.pop(i);
                end
            end
            
        end
        
        
        
        
        
        
        
        
        
        % Reset the Algorithm
        function this = reset(this)
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
                this.pop = this.sort_population(this.pop);
                this.best_sol = this.pop(1);
            end
            
        end
        
        % Sort Population
        function [pop, sort_order, obj_values] = sort_population(this, pop)
            
            % Check for Defined Optimization Problem
            direction = this.eff_problem.sort_direction;
            
            % Sort the Objective Values Vector
            [obj_values, sort_order] = sort([pop.obj_value], direction);
            
            % Sort (Re-order) Population
            pop = pop(sort_order);
            
        end
        
        % Sort and Select the Population
        function pop = sort_and_select(this, pop)
            
            % Sort Population
            pop = this.sort_population(pop);
            
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
            pos = reshape([pop.position], this.problem.var_count, [])';
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


function x = uniform_rand(lb, ub, varargin)
    % Generate Uniformly Distributed Random Numbers
    
    if ~exist('lb', 'var')
        lb = 0;
    end
    if ~exist('ub', 'var')
        ub = 1;
    end
    
    if isempty(varargin)
        mm = lb + ub;
        varargin{1} = size(mm);
    end
    
    x = lb + (ub - lb).*rand(varargin{:});
    
end