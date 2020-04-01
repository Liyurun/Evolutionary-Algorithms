function position =simulated_annealing(obj, problem)

opt = struct('max_iter', 10,...
                'pop_size' , 20,...
                'max_sub_iter',10,...
                'move_count', 15,...
                'mutation_rate', 0.0333,...
                'mutation_step_size', 0.5000,...
                'mutation_step_size_damp', 1,...
                'initial_temp', 100,...
                'final_temp',1,...
                'disp',1,...
                'problem',problem);
        sa = struct('name','SA');
        
        sa = opt;
        sa = start(sa,opt.max_iter);
            if opt.disp
                disp(['SA  started ...']);
                disp('Initializing population.');
            end
        sa = initialize(sa);
        
    
                % Iterations (Main Loop)
                for it = 1:sa.max_iter
                    
                    % Set Iteration Counter
                    sa.iter = it;
                    
                    % Iteration Started

                    sa = iterate(sa);
                    
                    % Update Histories
                    sa.nfe_history(it) = sa.nfe;
                    sa.best_obj_value_history(it) = sa.best_sol.obj_value;
                    position = sa.best_sol.position;
                    % Display Iteration Information
                    if sa.disp
                        disp(['Iteration ' num2str(sa.iter) ...
                              ': Best this. Value = ' num2str(sa.best_sol.obj_value) ...
                              ', position = ' num2str(position)]);
                    end
                    
 
                    
        
end
end



% Initialization
function sa = initialize(sa)

            % Create Initial Population (Sorted)
            sorted = true;
            sa  = init_pop(sa,sorted);
            
            % Initial Temperature
            sa.params.temp = sa.initial_temp;
            
            % Calculate Temperature Damp Rate
            sa.params.temp_damp = (sa.final_temp/sa.initial_temp)^(1/sa.max_iter);

            % Initial Value of Step Size
            sa.params.sigma = sa.mutation_step_size;
            
            
end


        

  % Initialize Population
        function sa = init_pop(sa, sorted)
            
            % Check for Sorted flag (Default is false, not sorted)
            if ~exist('sorted', 'var') || isempty(sorted)
                sorted = false;
            end
            
            % Initialize Population Array
            sa.pop = repmat(sa.empty_individual, sa.pop_size, 1);
            
            % Initialize Best Solution to the Worst Possible Solution
            sa.best_sol = sa.empty_individual;
            sa.best_sol.obj_value = sa.problem.worst_value;
            
            % Generate New Individuals
            for i = 1:sa.pop_size
                tpso = sa;
                % Generate New Solution
                sa.pop(i) = new_individual(tpso);
                
                % Compare to the Best Solution Ever Found
                if  sa.pop(i).obj_value < sa.best_sol.obj_value
                    sa.best_sol = sa.pop(i);
                end
                
            end
            

            
        end
% Iterations
function sa = iterate(sa)
            
            % Sub-Iterations
            for sub_iter = 1:sa.max_sub_iter
                
                % Create New Population
                newpop = repmat(sa.empty_individual, sa.pop_size, sa.move_count);
                for i = 1:sa.pop_size
                    for j = 1:sa.move_count
                        
                        % Perform Mutation (Move)
                        x = mutate(sa,sa.pop(i).position);
                        
                        % Evaluation
                        newpop(i,j) = new_individual(sa,x);
                        
                    end
                end
                
                % Columnize and Sort Newly Created Population
                newpop = sort_population(sa,newpop(:));
                
                % Compare the Best New Individual to Best Solution
                if newpop(1).obj_value < sa.best_sol.obj_value
                    sa.best_sol = newpop(1);
                end
                
                % Randomized Selection
                for i = 1:sa.pop_size
                    
                    % Check if new solution is better than current
                    if is_better(newpop(i), sa.pop(i))
                        
                        % If better, replace the old one
                        sa.pop(i) = newpop(i);
                        
                    else
                        
 
                            delta = newpop(i).obj_value - sa.pop(i).obj_value;
 
                        
                        % Compute Acceptance Probability
                        p = exp(-delta/sa.params.temp);
                        
                        % Accept / Reject
                        if rand() <= p
                            sa.pop(i) = newpop(i);
                        end
                        
                    end
                end
                
            end
            
            % Update Temperature
            sa.params.temp = sa.params.temp_damp * sa.params.temp;
            
            % Damp Step Size
            sa.params.sigma = sa.mutation_step_size_damp * sa.params.sigma;
            
        end
   
 %Perform Mutation (Create Neighbor Solution)
        function y = mutate(sa, x)
            
            % Mutation Rate
            mu = sa.mutation_rate;
            
            % Select Mutating Variables
            flag = (rand(size(x)) <= mu);
            if ~any(flag)
                % Select at least one variable to mutate
                j0 = randi(numel(x));
                flag(j0) = true;
            end
            j = find(flag);
            
            % Create Mutated Vector
            y = x;
            y(j) = x(j) + sa.params.sigma*randn(size(j));
            
            % Clip the Output
            y = ypea_clip(y, 0, 1);
            
        end
        
        
% Create a New Individual
        function ind = new_individual(sa, x)
            
            if ~exist('x', 'var') || isempty(x)
                x = rand(1,sa.problem.dim);
            end
            
            ind = sa.empty_individual;
            ind.position = clip(x, 0, 1);
            [ind.obj_value, ind.solution] = decode_and_eval(sa,x);
        end
        
        % Decode and Evaluate a Coded Solution
        function [z, sol] = decode_and_eval(sa, xhat)
        
        % Increment NFE
        sa.nfe = sa.nfe + 1;
        
        % Decode and Evaluate
        if 1
            [z, sol] = sa.problem.func(xhat);
        else
            z = 0;
            sol = [];
        end
        end
 % Sort Population
        function [pop, sort_order, obj_values] = sort_population(sa, pop)
            
 
            % Sort the Objective Values Vector
            [obj_values, sort_order] = sort([pop.obj_value] );
            
            % Sort (Re-order) Population
            pop = pop(sort_order);
            
        end
   % Check if a solution is better than other
        function b = is_better(sa, x1, x2)
            b = x1.obj_value < x1.obj_value;
        end
                      
     

function pop = eval(pso, pop)
for i = 1:numel(pop)
    pop(i).position = clip(pop(i).position, 0, 1);
    [pop(i).obj_value, pop(i).solution] =  decode_and_eval(pso,pop(i).position);
end
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


function opt = start(opt,max_iter)
          opt.pop=[];
          opt.best_sol=[];
          opt.iter = 0;
          opt.nfe = 0;
          opt.nfe_history = nan(1, max_iter);
          opt.best_obj_value_history = nan(1, max_iter);
          opt.run_time = 0;
          opt.avg_eval_time = 0;
          opt.must_stop = false;
          opt.best_obj_value = inf;
          % Initialize Emprt Prticle (Individual)
            opt.empty_individual = [];
            opt.empty_individual.position = [];
            opt.empty_individual.velocity = [];
            opt.empty_individual.obj_value = [];
            opt.empty_individual.solution = [];
            opt.empty_individual.best.position = [];
            opt.empty_individual.best.obj_value = [];
            opt.empty_individual.best.solution = [];
          
end



function opt = use_constriction_coeffs(opt)
            phi1 = opt.phi1
            phi2 = opt.phi2
            % Checing values of phi1 and phi2
            
            if ~exist('phi1', 'var')
                phi1 = [];
            else                
                validateattributes(phi1, {'numeric'}, {'scalar', 'nonnegative'});
            end
            
            if ~exist('phi2', 'var')
                phi2 = [];
            else
                validateattributes(phi2, {'numeric'}, {'scalar', 'nonnegative'});
            end
            
            if ~isempty(phi1)
                if ~isempty(phi2)
                    if phi1 + phi2 == 0
                        phi1 = 2.05;
                        phi2 = 2.05;
                    end
                else
                    phi2 = phi1;
                end
            else
                if ~isempty(phi2)
                    phi1 = phi2;
                else
                    phi1 = 2.05;
                    phi2 = 2.05;
                end
            end
            
            % Ensure than phi1 + phi2 is greater than or equal to 4.
            if phi1 + phi2 < 4
                phi = phi1 + phi2;
                phi1 = phi1/phi*4;
                phi2 = phi2/phi*4;
            end
            
            % Calculate Constriction Coefficients
            phi = phi1 + phi2;
            chi = 2 / (phi - 2 + sqrt(phi^2 - 4*phi));
            opt.w = chi;
            opt.wdamp = 1;
            opt.c1 = chi*phi1;
            opt.c2 = chi*phi2;
        end
