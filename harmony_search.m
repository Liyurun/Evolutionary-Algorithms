function position = harmony_search(opt, problem)


opt = struct('max_iter', 10,...
                'pop_size' , 20,...
                'hms'  ,20,...
                'new_count',  40,...
                'hmcr',  0.95,...
                'par',  0.1,...
                'fret_width',  0.1,...
                'fret_width_damp',  0.998,...
                'disp',1,...
                'problem',problem);
        hs = struct('name','BBO');
        hs = opt;
        hs = start(hs);
        
            if hs.disp
                disp(['BBO  ' 'started ...']);
                disp('Initializing population.');
            end       
        
        hs = initialize(hs);
        
        % Iterations (Main Loop)
                for it = 1:hs.max_iter
                    
                    % Set Iteration Counter
                    hs.iter = it;
                    
                    % Iteration Started
                    hs = iterate(hs);
                    
                    % Update Histories
                    hs.position_history(it).position = hs.best_sol.position;
                    hs.nfe_history(it) = hs.nfe;
                    hs.best_obj_value_history(it) = hs.best_sol.obj_value;
                    
          
                    % Display Iteration Information
                    if hs.disp
                        disp(['Iteration ' num2str(hs.iter) ...
                              ': Best de. Value = ' num2str(hs.best_sol.obj_value) ...
                              ', position = ' num2str(hs.best_sol.position)]);
                    end
                    
                    position =hs.best_sol.position;

                end
end


 % Initialization
        function this = initialize(this)
            sorted = true;
            this = init_pop(this,sorted);
            this.params.fw = this.fret_width;
        end
        
        % Iteration
        function this = iterate(this)
            
            % Decision Vector Size
            var_size = [1 this.problem.dim];
            
            % Decision Variable Count
            var_count = this.problem.dim;
            
            % Initialize Array for New Harmonies
            newpop = repmat(this.empty_individual, this.new_count, 1);
            
            % Create New Harmonies
            for k = 1:this.new_count

                % Create New Harmony Position
                newpop(k).position = rand(var_size);
                
                % Adjust New Harmony
                for j = 1:var_count
                    
                    % Use Harmony Memory
                    if rand <= this.hmcr
                        
                        i = randi([1 this.hms]);
                        newpop(k).position(j) = this.pop(i).position(j);
                        
                        % Pitch Adjustment
                        if rand <= this.par
                            delta = this.fret_width * randn();
                            newpop(k).position(j) = newpop(k).position(j) + delta;
                        end
                        
                    end
                    
                end
                
                % Clip Solution
                newpop(k).position = clip(newpop(k).position, 0, 1);
                
                % Evaluation
                newpop(k) = eval(this,newpop(k));
                
            end

            % Merge, Sort and Selection
            this.pop = sort_and_select(this,[this.pop; newpop]);
            
            % Determine Best Solution
            this.best_sol = this.pop(1);
            
            % Damp Fret Width
            this.params.fw = this.fret_width_damp * this.params.fw;
            
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
                [pop(i).obj_value, pop(i).solution] = decode_and_eval(this,pop(i).position);
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
