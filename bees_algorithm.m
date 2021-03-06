function position = bees_algorithm(obj, problem)


    opt = struct('max_iter', 10,...
                'pop_size' , 20,...
                'selected_site_ratio' , 0.5,...
                'selected_site_bee_ratio' , 0.1,...
                'selected_site_count', 15,...
                'selected_site_bee_count',3,...
                'scout_bee_count',30,...
                'elite_site_ratio' , 0.4,...
                'elite_site_bee_ratio' , 2,...
                'elite_site_bee_count',6,...
                'elite_site_count',6,...
                'dance_radius', 0.1,...
                'dance_radius_damp', 0.99,...
                'disp',1,...
                'problem',problem);
        ba = struct('name','BA');
        ba = opt;
        ba = start(ba);
        
            if ba.disp
                disp(['BA  ' 'started ...']);
                disp('Initializing population.');
            end       
        
        ba = initialize(ba);
        
        % Iterations (Main Loop)
                for it = 1:ba.max_iter
                    
                    % Set Iteration Counter
                    ba.iter = it;
                    
                    % Iteration Started
                    ba = iterate(ba);
                    
                    % Update Histories
                    ba.position_history(it).position = ba.best_sol.position;
                    ba.nfe_history(it) = ba.nfe;
                    ba.best_obj_value_history(it) = ba.best_sol.obj_value;
                    
          
                    % Display Iteration Information
                    if ba.disp
                        disp(['Iteration ' num2str(ba.iter) ...
                              ': Best de. Value = ' num2str(ba.best_sol.obj_value) ...
                              ', position = ' num2str(ba.best_sol.position)]);
                    end
                    
                    position =ba.best_sol.position;

                end
                

        
end

        % Initialization
        function this = initialize(this)
            
            % Create Initial Population (Sorted)
            sorted = true;
            this = init_pop(this,sorted);
            
            % Initial Value of Dance Radius
            this.params.r = this.dance_radius;
            
        end
        
        % Iterations
        function this = iterate(this)
            
            % Iterate based on algorithm type
            % Standard BA
                this = iterate_standard(this);
            
            % Sort Population
            this.pop = sort_population(this,this.pop);
            
            % Update Best Solution Ever Found
            this.best_sol = this.pop(1);
            
            % Damp Dance Radius
            this.params.r = this.dance_radius_damp * this.params.r;
            
        end
        
        % Standard BA Iterator
        function this = iterate_standard(this)
            
            % Elite Sites
            for i = 1:this.elite_site_count
                
                % Create New Bees (Solutions)
                best_new_bee.obj_value = inf;
                for j = 1:this.elite_site_bee_count
                    xnew = perform_dance(this,this.pop(i).position);
                    new_bee = new_individual(this,xnew);
                    if  is_better(this,new_bee, best_new_bee)
                        best_new_bee = new_bee;
                    end
                end
                
                % Compare to Best Solution Ever Found
                if is_better(this,best_new_bee, this.pop(i))
                    this.pop(i) = best_new_bee;
                end
                
            end

            % Selected Non-Elite Sites
            for i = this.elite_site_count+1:this.selected_site_count
                
                % Create New Bees (Solutions)
                best_new_bee.obj_value = inf;
                for j = 1:this.selected_site_bee_count
                    xnew = perform_dance(this,this.pop(i).position);
                    new_bee = new_individual(this,xnew);
                    if is_better(this,new_bee, best_new_bee)
                        best_new_bee = new_bee;
                    end
                end
                
                % Compare to Best Solution Ever Found
                if is_better(this,best_new_bee, this.pop(i))
                    this.pop(i) = best_new_bee;
                end
                
            end

            % Non-Selected Sites
            for i = this.selected_site_count+1:this.scout_bee_count
                this.pop(i) =  new_individual(this);
            end
            
        end
        
        % Perform Bee Dance
        function y = perform_dance(this, x)
            j = randi([1 numel(x)]);
            y = x;
            y(j) = x(j) + this.params.r * uniform_rand(-1, 1);
            y = clip(y, 0, 1);
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
