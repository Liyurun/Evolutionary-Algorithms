function position = genetic_algorithm(obj, problem)

        opt = struct('max_iter',10,...
                    'pop_size' , 20,...
                    'crossover_prob' , 0.7,...
                    'crossover_inflation' , 0.4,...
                    'mutation_prob' , 0.3,...
                    'mutation_rate' , 0.1,...
                    'mutation_step_size' , 0.5,...
                    'mutation_step_size_damp' , 0.99,...
                    'selection_pressure',5,...
                    'sigma',0.5,...
                    'disp',1,...
                    'problem',problem);
        ga = struct('name','GA');
        ga = opt;
        ga = start(ga);
        
            if ga.disp
                disp(['DE  ' 'started ...']);
                disp('Initializing population.');
            end       
        
        ga = initialize(ga);
        
                % Iterations (Main Loop)
                for it = 1:ga.max_iter
                    
                    % Set Iteration Counter
                    ga.iter = it;
                    % Iteration Started
                    ga = iterate(ga);
                    
                    % Update Histories
                    ga.nfe_history(it) = ga.nfe;
                    ga.best_obj_value_history(it) = ga.best_sol.obj_value;
                    
 
                    
                    % Display Iteration Information
                    if ga.disp
                        disp(['Iteration ' num2str(ga.iter) ...
                              ': Best this. Value = ' num2str(ga.best_sol.obj_value) ...
                              ', position = ' num2str([ga.best_sol.position])]);
                    end
                    
                    position = ga.best_sol.position;

                end
                
end
function ga = iterate(ga)
 
            % Population Size
            pop_size = ga.pop_size;

                beta = ga.selection_pressure;
                obj_values = ga.pop.obj_value;
                P =  get_selection_probs(ga,obj_values, beta);

            % Perform Crossover
            nc = 2*round(ga.crossover_prob * pop_size/2);
            popc = repmat(ga.empty_individual, nc/2, 2);
            for k = 1:nc/2
                
                % Select Parents
                    i1 = roulette_wheel_selection(P);
                    i2 = roulette_wheel_selection(P);
                % Perform Crossover
                x1 = ga.pop(i1).position;
                x2 = ga.pop(i2).position;
                [y1, y2]=crossover(ga,x1, x2);
                
                % Evaluate Offsprings
                popc(k,1) = new_individual(ga,y1);
                popc(k,2) = new_individual(ga,y2);
                
            end
            popc = popc(:);
            
            % Perform Mutation
            nm = round(ga.mutation_prob * pop_size);
            popm = repmat(ga.empty_individual, nm, 1);
            for k = 1:nm
                
                % Select Parent
                i = randi([1 pop_size]);
                x = ga.pop(i).position;
                
                % Perform Mutation
                y =  mutate(ga,x);
                
                % Evaluate Offspring
                popm(k) = new_individual(ga,y);
                
            end
            
            % Merge, Sort and Selection
            ga.pop = sort_and_select(ga,[ga.pop; popm; popc]);
            
            % Update Best Solution Ever Found
            ga.best_sol = ga.pop(1);
            
            % Damp Mutation Step Size
            ga.sigma = ga.mutation_step_size_damp * ga.sigma;          
end
function y = mutate(ga, x)
            
            mu = ga.mutation_rate;
            sigma = ga.sigma;
            
            nvar = numel(x);
            nmu = ceil(mu * nvar);
            
            j = rand_sample(nvar, nmu);
            
            y = x;
            y(j) = x(j) + sigma*randn(size(j));
            
            y = clip(y, 0, 1);
            
        end
function [y1, y2] = crossover(ga, x1, x2)
            
            gamma = ga.crossover_inflation;
            alpha = uniform_rand(-gamma, 1+gamma, size(x1));
            
            y1 = alpha.*x1 + (1-alpha).*x2;
            y2 = alpha.*x2 + (1-alpha).*x1;
            
            y1 = clip(y1, 0, 1);
            y2 = clip(y2, 0, 1);
            
        end
function p = get_selection_probs(ga, values, selection_pressure)
            
            % Check if Selection Pressure specified
            if ~exist('selection_pressure', 'var')
                selection_pressure = 1;
            end
            
            % Change selection pressure sign, according to problem type
 
                alpha = -selection_pressure;
 
            
            % Calculate Selection Probabilities
            p = exp(alpha*values);
            p = p/sum(p);
            
        end
function pop = sort_and_select(ga, pop)
            
            % Sort Population

              [obj_values, sort_order] = sort([pop.obj_value]);
              pop = pop(sort_order);
            % Set the Population Size Limit
            n = min(ga.pop_size, numel(pop));
            pop = pop(1:n);
            
        end
                
    function ga = start(ga)
        ga.empty_individual.position = [];
        ga.empty_individual.obj_value = [];
        ga.empty_individual.solution = [];
        
        ga.pop= [];
        ga.best_sol = [];
        ga.iter = 0;
        ga.nfe = 0;
        ga.run_time = 0;
        ga.avg_eval_time = 0;
        ga.nfe_history = [];
        ga.best_obj_value_history = [];
        ga.iter = 0;
        % specitial params
        ga.nfe = 0;
        ga.nfe_history = nan(1, ga.max_iter);
        ga.best_obj_value_history = nan(1, ga.max_iter);
        ga.run_time = 0;
        ga.avg_eval_time = 0;
        ga.must_stop = false;
    end
    
function ga = initialize(ga)
sorted = false;

% Check for Sorted flag (Default is false, not sorted)
if ~exist('sorted', 'var') || isempty(sorted)
    sorted = false;
end

% Initialize Population Array
ga.pop = repmat(ga.empty_individual, ga.pop_size, 1);

% Initialize Best Solution to the Worst Possible Solution
ga.best_sol = ga.empty_individual;
ga.best_sol.obj_value = inf;

% Generate New Individuals
for i = 1:ga.pop_size
    
    % Generate New Solution
    x = rand(1,ga.problem.dim);
    ga.pop(i) = new_individual(ga,x);
    
    % Compare to the Best Solution Ever Found
    if ga.pop(i).obj_value< ga.best_sol.obj_value
        ga.best_sol = ga.pop(i);
    end
    
end

% Sort the Population if it is needed
if sorted
    ga.pop = ga.sort_population(ga.pop);
    ga.best_sol = ga.pop(1);
end

end
function ind = new_individual(de,x)

ind = de.empty_individual;
ind.position = x;
[ind.obj_value, ind.solution] = decode_and_eval(de,x);
end
function [z, sol] = decode_and_eval(de, x)

% Increment NFE
de.nfe = de.nfe + 1;

% Decode and Evaluate
[z, sol] = de.problem.func(x);
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