function position = ant_colony_optimization(obj, problem)

    opt = struct('max_iter', 10,...
                'pop_size' , 20,...
                'sample_count',50,...
                'q',0.5,...
                'zeta',1,...
                'disp',1,...
                'problem',problem);
        aco = struct('name','ACO');
        aco = opt;
        aco = start(aco);
        
            if aco.disp
                disp(['ACO  ' 'started ...']);
                disp('Initializing population.');
            end       
        
        aco = initialize(aco);
        
        % Iterations (Main Loop)
                for it = 1:aco.max_iter
                    
                    % Set Iteration Counter
                    aco.iter = it;
                    
                    % Iteration Started
                    aco = iterate(aco);
                    
                    % Update Histories
                    aco.position_history(it).position = aco.best_sol.position;
                    aco.nfe_history(it) = aco.nfe;
                    aco.best_obj_value_history(it) = aco.best_sol.obj_value;
                    
          
                    % Display Iteration Information
                    if aco.disp
                        disp(['Iteration ' num2str(aco.iter) ...
                              ': Best de. Value = ' num2str(aco.best_sol.obj_value) ...
                              ', position = ' num2str(aco.best_sol.position)]);
                    end
                    
                    position =aco.best_sol.position;

                end
                

        
        
end




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

 % Initialization
 function this =  initialize(this)

            % Create Initial Population (Sorted)
            sorted = true;
            this =init_pop(this,sorted);
            
            % Calculate Selection Probabilities
            qn = this.q * this.pop_size;
            w = 1/(sqrt(2*pi)*qn)*exp(-0.5*(((1:this.pop_size)-1)/qn).^2);
            this.params.p = normalize_probs(w);
            
        end
        
        % Iterations
        function this = iterate(this)
            
            % Decision Vector Size
            var_size = [1 this.problem.dim];
            
            % Means and Standard Deviations
            s = get_positions(this,this.pop);
            sigma = zeros(size(s));
            for l = 1:this.pop_size
                D = sum(abs(repmat(s(l,:), this.pop_size, 1) - s));
                sigma(l,:) = this.zeta*D/(this.pop_size - 1);
            end
            
            % Generate Samples
            new_pop = repmat(this.empty_individual, this.sample_count, 1);
            for k = 1:this.sample_count
                x = zeros(var_size);
                for i = 1:numel(x)
                    l = roulette_wheel_selection(this.params.p);
                    x(i) = s(l,i) + sigma(l,i)*randn();
                end
                new_pop(k) = new_individual(this,x);
            end
            
            % Merge, Sort and Selection
            this.pop = sort_and_select(this,[this.pop; new_pop]);
            
            % Determine Best Solution
            this.best_sol = this.pop(1);
            
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






