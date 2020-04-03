function position = future_search_algorithm(obj, problem)

        opt = struct('max_iter', 10,...
                        'pop_size' , 20,... 
                        'r_time',10,...	 
                        'disp',1,...
                        'problem',problem);
        fsa = struct('name','FSA');
        fsa = opt;
        fsa = start(fsa);
        
            if fsa.disp
                disp(['FSA  ' 'started ...']);
                disp('Initializing population.');
            end       
        
        fsa = initialize(fsa);
        
        % Iterations (Main Loop)
                for it = 1:fsa.max_iter
                    
                    % Set Iteration Counter
                    fsa.iter = it;
                    
                    % Iteration Started
                    fsa = iterate(fsa);
                    
                    % Update Histories
                    fsa.position_history(it).position = fsa.best_sol.position;
                    fsa.nfe_history(it) = fsa.nfe;
                    fsa.best_obj_value_history(it) = fsa.best_sol.obj_value;
                    
          
                    % Display Iteration Information
                    if fsa.disp
                        disp(['Iteration ' num2str(fsa.iter) ...
                              ': Best de. Value = ' num2str(fsa.best_sol.obj_value) ...
                              ', position = ' num2str(fsa.best_sol.position)]);
                    end
                    
                    position =fsa.best_sol.position;

                end

end



        % Initialization
        function this = initialize(this)
            
            % Create Initial Population (Not Sorted)
            sorted = false;
            this = init_pop(this,sorted);
            for ii = 1: this.pop_size
                this.Lbe(ii).obj_value = this.pop(ii).obj_value
                this.Lbep(ii).position = this.pop(ii).position
            end
            % Initial Value of Mutation Coefficient
            %this.params.alpha = this.alpha;
            
        end
        
        % Iterations
        function this = iterate(this)
            
          
            % Decision Vector Size
            var_size = [1 this.problem.dim];
            % Create New Population
            newpop = repmat(this.empty_individual, this.pop_size, 1);
            for i = 1:this.pop_size
                
                % Initialize to Worst Objective Value
                newpop(i).obj_value = this.problem.worst_value;
                % Calculate objective function for each search agent

                this.pop(i).position=this.pop(i).position+(-this.pop(i).position+this.best_sol.position)*rand+(-this.pop(i).position+this.Lbep(i).position)*rand;
                for k=1:this.problem.dim
                    if this.pop(i).position(k)>this.problem.upbound(k), this.pop(i).position(k)=this.problem.upbound(k); end
                    if this.pop(i).position(k)<this.problem.lobound(k), this.pop(i).position(k)=this.problem.lobound(k); end
                end               
                this.pop(i).obj_value=this.problem.func(this.pop(i).position);
                
                % Update the loacal best solution
                       if (this.pop(i).obj_value<=this.Lbe(i).obj_value) 
                            this.Lbep(i).position=this.pop(i).position;
                            this.Lbe(i).obj_value=this.pop(i).obj_value;
                       end
                
                if this.pop(i).obj_value<this.best_sol.obj_value % Change this to > for maximization problem
                    this.best_sol.obj_value=this.pop(i).obj_value; % Update alpha
                    this.best_sol.position=this.pop(i).position;
                end
                
            end
            

            % loop of  the initial update
            for i=1:this.pop_size,
                Si(i,:) = this.best_sol.position + (this.best_sol.position - this.Lbep(i).position).*rand;
                 
                for k=1:this.problem.dim
                    if Si(i,k)>this.problem.upbound(k), Si(i,k)=this.problem.upbound(k); end
                    if Si(i,k)<this.problem.lobound(k), Si(i,k)=this.problem.lobound(k); end
                end  
                 this.pop(i).obj_value=this.problem.func(Si(i,:));
            % Update the loacal best solution
                 if ( this.pop(i).obj_value<=this.pop(i).obj_value) 
                    this.pop(i).position=Si(i,:);
                    this.Lbep(i).position=Si(i);
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
        
 