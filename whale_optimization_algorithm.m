function position = whale_optimization_algorithm(obj, problem)


opt = struct('max_iter', 100,...
                'pop_size' , 200,...
                'SearchAgents_no',30,...
                'disp',1,...
                'problem',problem);
        woa = struct('name','WOA');
        woa = opt;
        woa = start(woa);
        
            if woa.disp
                disp(['WOA  ' 'started ...']);
                disp('Initializing population.');
            end       
        
        woa = initialize(woa);
        
        % Iterations (Main Loop)
                for it = 1:woa.max_iter
                    
                    % Set Iteration Counter
                    woa.iter = it;
                    
                    % Iteration Started
                    woa = iterate(woa);
                    
                    % Update Histories
                    woa.position_history(it).position = woa.best_sol.position;
                    woa.nfe_history(it) = woa.nfe;
                    woa.best_obj_value_history(it) = woa.best_sol.obj_value;
                    
          
                    % Display Iteration Information
                    if woa.disp
                        disp(['Iteration ' num2str(woa.iter) ...
                              ': Best de. Value = ' num2str(woa.best_sol.obj_value) ...
                              ', position = ' num2str(woa.best_sol.position)]);
                    end
                    
                    position =woa.best_sol.position;

                end
                
end







        % Initialization
        function this = initialize(this)
            
            % Create Initial Population (Not Sorted)
            sorted = false;
            this = init_pop(this,sorted);
            
            % Initial Value of Mutation Coefficient
            %this.params.alpha = this.alpha;
            
        end
        
        % Iterations
        function this = iterate(this)
            
          
            % Decision Vector Size
            var_size = [1 this.problem.dim];
            this.best_sol.obj_value = inf;
            this.best_sol.position = [];
            % Create New Population
            newpop = repmat(this.empty_individual, this.pop_size, 1);
            for i = 1:this.pop_size
                
                % Initialize to Worst Objective Value
                newpop(i).obj_value = this.problem.worst_value;
                % Calculate objective function for each search agent
                    
                this.pop(i).obj_value=this.problem.func(this.pop(i).position);
                if this.pop(i).obj_value<this.best_sol.obj_value % Change this to > for maximization problem
                    this.best_sol.obj_value=this.pop(i).obj_value; % Update alpha
                    this.best_sol.position=this.pop(i).position;
                end
                
            end
            
            
    a=2-this.iter *((2)/this.max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+this.iter*((-1)/this.max_iter);
    
    % Update the Position of search agents 
    for i=1:size(this.pop,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(this.pop(1).position,2)
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(this.SearchAgents_no*rand()+1);
                    X_rand = this.pop(rand_leader_index).position(:);
                    D_X_rand=abs(C*X_rand(j)-this.pop(i).position(j)); % Eq. (2.7)
                    this.pop(i).position(j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    
                elseif abs(A)<1
                    D_Leader=abs(C*this.best_sol.position(j)-this.pop(i).position(j)); % Eq. (2.1)
                    this.pop(i).position(j)=this.best_sol.position(j)-A*D_Leader;      % Eq. (2.2)
                end
                
            elseif p>=0.5
              
                distance2Leader=abs(this.best_sol.position(j)-this.pop(i).position(j));
                % Eq. (2.5)
                this.pop(i).position(j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+this.best_sol.position(j);
                
            end
            
        end
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






