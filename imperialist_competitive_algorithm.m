function position = imperialist_competitive_algorithm(obj, problem)


opt = struct('max_iter', 100,...
                'pop_size' , 200,...
                'empire_count' , 5,...
                'selection_pressure' , 1,...
                'assimilation_coeff' , 1.5,...
                'revolution_prob' , 0.05,...
                'revolution_rate' , 0.1,...
                'revolution_step_size' , 0.1,...
                'revolution_step_size_damp' , 0.99,...
                'zeta' , 0.2,...
                'colony_count',90,...
                'disp',1,...
                'problem',problem);
        ica = struct('name','ICA');
        ica = opt;
        ica = start(ica);
        
            if ica.disp
                disp(['ICA  ' 'started ...']);
                disp('Initializing population.');
            end       
        
        ica = initialize(ica);
        
        % Iterations (Main Loop)
                for it = 1:ica.max_iter
                    
                    % Set Iteration Counter
                    ica.iter = it;
                    
                    % Iteration Started
                    ica = iterate(ica);
                    
                    % Update Histories
                    ica.position_history(it).position = ica.best_sol.position;
                    ica.nfe_history(it) = ica.nfe;
                    ica.best_obj_value_history(it) = ica.best_sol.obj_value;
                    
          
                    % Display Iteration Information
                    if ica.disp
                        disp(['Iteration ' num2str(ica.iter) ...
                              ': Best de. Value = ' num2str(ica.best_sol.obj_value) ...
                              ', position = ' num2str(ica.best_sol.position)]);
                    end
                    
                    position =ica.best_sol.position;

                end


end




        % Initialization
        function this = initialize(this)
            
            % Create Initial Population (Sorted)
            sorted = true;
            this = init_pop(this,sorted);
            
            % Determine Imperialists and Colonies
            pop = this.pop;
            imp = pop(1:this.empire_count);
            col = pop(this.empire_count+1:end);
            
            % Empty Empire
            empty_empire.imp = [];
            empty_empire.col = repmat(this.empty_individual, 0, 1);
            empty_empire.colony_count = 0;
            empty_empire.total_obj_value = [];
            
            % Initialize Empires
            this.emp = repmat(empty_empire, this.empire_count, 1);
            
            % Assign Imperialists
            for k = 1:this.empire_count
                this.emp(k).imp = imp(k);
            end
            
            % Determine Selection Probabilities
            obj_values = [imp.obj_value];
            obj_values = obj_values/max(abs(obj_values));
            P = get_selection_probs(this,obj_values, this.selection_pressure);
            
            % Assign Colonies
            for j = 1:this.colony_count
                k = roulette_wheel_selection(P);
                this.emp(k).col = [this.emp(k).col; col(j)];
                this.emp(k).colony_count = numel(this.emp(k).col);
            end
            
            % Initial Value of Step Size
            this.params.sigma = this.revolution_step_size;
            
            % Update Total Objective Values of Empires
            this = update_empires_total_objective_values(this);
            
        end
        
        function this = iterate(this)
            
            % Assimilation
            this = assimilation(this);
            
            % Revolution
            this = revolution(this);
            
            % Intra-Empire Competition
            this = intra_empire_competition(this);

            % Update Total Objective Values of Empires
            this = update_empires_total_objective_values(this);

            % Inter-Empire Competition
            this = inter_empire_competition(this);
            
            % Update Revolution Step Size
            this.params.sigma = this.revolution_step_size_damp * this.params.sigma;
            
        end
        
        function this = update_empires_total_objective_values(this)

            for k = 1:numel(this.emp)
                if this.emp(k).colony_count > 0
                    this.emp(k).total_obj_value = this.emp(k).imp.obj_value ...
                        + this.zeta*mean([this.emp(k).col.obj_value]);
                else
                    this.emp(k).total_obj_value = this.emp(k).imp.obj_value;
                end
            end

        end        
        
        function this = assimilation(this)
            
            % Decision Vector Size
            var_size = [1 this.problem.dim];
            
            % Assimilation Coefficient
            beta = this.assimilation_coeff;
            
            % Assimilate Colonies
            for k = 1:numel(this.emp)
                for i = 1:this.emp(k).colony_count
                    
                    % Create New Solution
                    this.emp(k).col(i).position = this.emp(k).col(i).position ...
                        + beta*rand(var_size).*(this.emp(k).imp.position - this.emp(k).col(i).position);
                    
                    % Evaluation
                    this.emp(k).col(i) = eval(this,this.emp(k).col(i));
                    
                    % Compare to Best Solution Ever Found
                    if is_better(this, this.emp(k).col(i), this.best_sol)
                        this.best_sol = this.emp(k).col(i);
                    end
                    
                end
            end
            
        end
        
        function this = revolution(this)
            
            % Decision Vector Size
            var_size = [1 this.problem.dim];
            
            % Decision Variables Count
            var_count = this.problem.dim;
            
            % Revolution Rate
            mu = this.revolution_rate;
            
            % Number of Revoluted Decision Variables
            nmu = ceil(mu*var_count);
            
            % Revolution Step Size
            sigma = this.params.sigma;
            
            % Revolve Imperialists and Colonies
            for k = 1:numel(this.emp)
                
                % Apply Revolution to Imperialist
                xnew = this.emp(k).imp.position + sigma*randn(var_size);
                jj = rand_sample(var_count, nmu);
                newimp = this.emp(k).imp;
                newimp.position(jj) = xnew(jj);
                newimp = eval(this,newimp);
                if is_better(this, newimp, this.emp(k).imp)
                    this.emp(k).imp = newimp;
                    if is_better(this, this.emp(k).imp, this.best_sol)
                        this.best_sol = this.emp(k).imp;
                    end
                end
                
                % Apply Revolution to Colonies
                for i = 1:this.emp(k).colony_count
                    if rand <= this.revolution_prob
                        xnew = this.emp(k).col(i).position + sigma*randn(var_size);
                        jj = rand_sample(var_count, nmu);
                        this.emp(k).col(i).position(jj) = xnew(jj);
                        this.emp(k).col(i) = eval(this,this.emp(k).col(i));
                        if is_better(this, this.emp(k).col(i), this.best_sol)
                            this.best_sol = this.emp(k).col(i);
                        end
                    end
                end
            end
            
        end
        
        function this = intra_empire_competition(this)
            
            for k = 1:numel(this.emp)
                for i = 1:this.emp(k).colony_count
                    
                    % Compare Colonies of Empires to Corresponding Imperialist
                    if is_better(this, this.emp(k).col(i), this.emp(k).imp)
                        
                        % If colony is better, then swap colony and imp.
                        
                        imp = this.emp(k).imp;
                        col = this.emp(k).col(i);
                        
                        this.emp(k).imp = col;
                        this.emp(k).col(i) = imp;
                        
                    end
                    
                end
            end
            
        end

        function this = inter_empire_competition(this)

            % In case of one empire, inter-empire competition is not needed
            if numel(this.emp) == 1
                return;
            end
            
            % Gather and Normalize Total Objective Values
            total_obj_values = [this.emp.total_obj_value];
            total_obj_values = total_obj_values/max(abs(total_obj_values));
            
            % Find Weakest Empire
            [~, weakest_empire_index] = min(total_obj_values);
            weakest_empire = this.emp(weakest_empire_index);
            
            % Calculate Selection Probabilities
            P = get_selection_probs(this,total_obj_values, this.selection_pressure);
            P(weakest_empire_index)=0;
            P = normalize_probs(P);
            
            % If the weakset empire has any colonies
            if weakest_empire.colony_count > 0
                
                % Find Weakest Colony of th Weakest Empire
                [~, weakest_colony_index] = min([weakest_empire.col.obj_value]);
                weakest_col = weakest_empire.col(weakest_colony_index);
                
                winner_empire_index = roulette_wheel_selection(P);
                winner_empire = this.emp(winner_empire_index);

                winner_empire.col(end+1) = weakest_col;
                winner_empire.colony_count = numel(winner_empire.col);
                this.emp(winner_empire_index) = winner_empire;
                
                weakest_empire.col(weakest_colony_index)=[];
                weakest_empire.colony_count = numel(weakest_empire.col);
                this.emp(weakest_empire_index) = weakest_empire;
                
            end
            
            % If the weakest empire has no colonies
            if weakest_empire.colony_count == 0
                
                winner_empire_index_2 = roulette_wheel_selection(P);
                winner_empire_2 = this.emp(winner_empire_index_2);

                winner_empire_2.col(end+1) = weakest_empire.imp;
                this.emp(winner_empire_index_2) = winner_empire_2;

                this.emp(weakest_empire_index)=[];
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

        
        
        function s = rand_sample(n, k)
    % Randomly selecting k samples from n items
    
    if ~exist('k', 'var')
        k = 1;
    end

    k = min(n, k);
    
    r = rand(1,n);
    [~, so] = sort(r);
    
    s = so(1:k);
    
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
        
