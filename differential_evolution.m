function position = Differential_Evolution(obj, problem)

    opt = struct('max_iter', 10,...
                'pop_size' , 20,...
                'beta_min',0.1,...
                'beta_max',0.9,...
                'crossover_prob',0.1,...
                'disp',1,...
                'problem',problem);
        de = struct('name','DE');
        de = opt;
        de = start(de);
        
            if de.disp
                disp(['DE  ' 'started ...']);
                disp('Initializing population.');
            end       
        
        de = initialize(de);
        
        % Iterations (Main Loop)
                for it = 1:de.max_iter
                    
                    % Set Iteration Counter
                    de.iter = it;
                    
                    % Iteration Started
                    de = iterate(de);
                    
                    % Update Histories
                    de.position_history(it).position = de.best_sol.position;
                    de.nfe_history(it) = de.nfe;
                    de.best_obj_value_history(it) = de.best_sol.obj_value;
                    
          
                    % Display Iteration Information
                    if de.disp
                        disp(['Iteration ' num2str(de.iter) ...
                              ': Best de. Value = ' num2str(de.best_sol.obj_value) ...
                              ', position = ' num2str(de.best_sol.position)]);
                    end
                    
                    position =de.best_sol.position;

                end
                

        
        
        
end

function de = start(de)

de.base_vector= 'best';
de.diff_vectors_count= 2;
de.crossover_method= 'exp';
de.crossover_prob= 0.1000;
de.type= 'DE/best/2/exp';
de.empty_individual.position = [];
de.empty_individual.obj_value = [];
de.empty_individual.solution = [];
de.pop= [];
de.best_sol= [];
de.iter= 0;
de.nfe= 0;
de.run_time= 0;
de.avg_eval_time= 0;
de.nfe_history= [];
de.best_obj_value_history= [];
de.must_stop= 0;
% reset

            de.iter = 0;
            de.nfe = 0;
            de.nfe_history = nan(1, de.max_iter);
            de.best_obj_value_history = nan(1, de.max_iter);
            de.run_time = 0;
            de.avg_eval_time = 0;
            de.must_stop = false;

end

% Iterations
function de = iterate(de)

for i = 1:de.pop_size
    
    % Position of Selected Solution
    x = de.pop(i).position;
    
    % Perform Mutation
    y = perform_mutation(de,i);
    
    % Crossover
    z = perform_crossover(de,x, y);
    
    % Create and Evaluate New Solution
    newsol =  new_individual(de,z);
    
    % Compare to Current Solution
    if newsol.obj_value < de.pop(i).obj_value
        
        % Replace Current Solution
        de.pop(i) = newsol;
        
        % Compare to Best Solution ever found
        if de.pop(i).obj_value <  de.best_sol.obj_value
            de.best_sol = de.pop(i);
        end
    end
    
end

end
function y = perform_mutation(de, i)

% Decision Vector Size
var_size = de.problem.dim;

% Number of Difference Vectors
ndv = de.diff_vectors_count;

% Acceleration Coefficient;
beta =   de.beta_min + (de.beta_max - de.beta_min).*rand(1,var_size);



% Number of Individuals Needed

ndv = min(ndv, floor((de.pop_size - 1)/2));
ns = 2*ndv;

% Select Inidividuals
A = randperm(de.pop_size);
A(A==i) = [];
A = A(1:ns);

% Calculate Base Vector
xbase = de.best_sol.position;
% Calculate Differences and Final Mutated Vector
A = reshape(A, 2, []);
a = A(1,:);
b = A(2,:);
y = xbase;
for k = 1:numel(a)
    y = y + beta().*(de.pop(a(k)).position - de.pop(b(k)).position);
end

end

function z = perform_crossover(de, x, y)

n = numel(x);

switch de.crossover_method
    case 'bin'
        % Binomial Crossover
        J = (rand(size(x)) <= de.crossover_prob);
        J(randi(n)) = true;
        
    case 'exp'
        % Exponential Crossover
        S = randi(n);
        L = randi(n);
        J = mod(S:S+L, n);
        J(J==0) = n;
end

z = x;
z(J) = y(J);

end


function de = initialize(de)
sorted = false;

% Check for Sorted flag (Default is false, not sorted)
if ~exist('sorted', 'var') || isempty(sorted)
    sorted = false;
end

% Initialize Population Array
de.pop = repmat(de.empty_individual, de.pop_size, 1);

% Initialize Best Solution to the Worst Possible Solution
de.best_sol = de.empty_individual;
de.best_sol.obj_value = inf;

% Generate New Individuals
for i = 1:de.pop_size
    
    % Generate New Solution
    x = rand(1,de.problem.dim);
    de.pop(i) = new_individual(de,x);
    
    % Compare to the Best Solution Ever Found
    if de.pop(i).obj_value< de.best_sol.obj_value
        de.best_sol = de.pop(i);
    end
    
end

% Sort the Population if it is needed
if sorted
    de.pop = de.sort_population(de.pop);
    de.best_sol = de.pop(1);
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


        