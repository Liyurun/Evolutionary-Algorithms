function position = particle_swarm_optimization(obj, problem)
   opt = struct('max_iter', 10,...
                'pop_size' , 20,...
                'w',0.5,...
                'wdamp',1,...
                'c1',1,...
                'c2' , 2,...
                'phi1' , 2.05,...
                'phi2' , 2.05,...
                'disp',1,...
                'problem',problem)
        pso = struct('name','PSO');
        % Make Use of Constriction Coefficients
        opt = use_constriction_coeffs(opt);
        % Define parameters before calculation
        pso = opt;
        pso = start(pso,opt.max_iter);
            if opt.disp
                disp(['PSO  started ...']);
                disp('Initializing population.');
            end
        pso = initialize(pso);
            
     
for it = 1:pso.max_iter
    
    % Set Iteration Counter
    pso.iter = it;
    
    % Iteration Started
   
    pso = iterate(pso);
    
    % Update Histories
    pso.position_history(it).position = pso.best_obj_value.position;
    pso.best_obj_value_history(it) = pso.best_obj_value.obj_value;
    
    % Iteration Ended
    
    
    % Display Iteration Information
    if pso.disp
        disp(['Iteration ' num2str(pso.iter) ...
            ': Best this. Value = ' num2str(pso.best_obj_value.obj_value) ...
            ', Best solution = ' num2str(pso.best_obj_value.position)]);
    end
    position = pso.best_obj_value.position;
    
%     % Check for Goals and Termination Conditions
%     pso.check_for_goals();
%     
%     % Check if it is needed to stop the execution
%     if pso.must_stop
%         break;
%     end
    
end
        
        
end


% Initialization
function pso = initialize(pso)

% Initialize Velocity Vector
pso.empty_individual.velocity = zeros(1,pso.problem.dim);

% Create Initial Population (Not Sorted)
sorted = false;
pso = init_pop(pso,sorted);

% Initialize Personal Bests
for i = 1:pso.pop_size
    pso.pop(i).best.position = pso.pop(i).position;
    pso.pop(i).best.solution = pso.pop(i).solution;
    pso.pop(i).best.obj_value = pso.pop(i).obj_value;
end

% Initial Value of Inertia Weight
pso.params.w = pso.w;

end

  % Initialize Population
        function pso = init_pop(pso, sorted)
            
            % Check for Sorted flag (Default is false, not sorted)
            if ~exist('sorted', 'var') || isempty(sorted)
                sorted = false;
            end
            
            % Initialize Population Array
            pso.pop = repmat(pso.empty_individual, pso.pop_size, 1);
            
            % Initialize Best Solution to the Worst Possible Solution
            pso.best_sol = pso.empty_individual;
            pso.best_sol.obj_value = pso.problem.worst_value;
            
            % Generate New Individuals
            for i = 1:pso.pop_size
                tpso = pso;
                % Generate New Solution
                pso.pop(i) = new_individual(tpso);
                
                % Compare to the Best Solution Ever Found
                if  pso.pop(i).obj_value < pso.best_sol.obj_value
                    pso.best_sol = pso.pop(i);
                end
                
            end
            

            
        end
        

        
% Create a New Individual
        function ind = new_individual(this, x)
            
            if ~exist('x', 'var') || isempty(x)
                x = rand(1,this.problem.dim);
            end
            
            ind = this.empty_individual;
            ind.position = clip(x, 0, 1);
            [ind.obj_value, ind.solution] = decode_and_eval(this,x);
        end
        
        % Decode and Evaluate a Coded Solution
        function [z, sol] = decode_and_eval(this, xhat)
        
        % Increment NFE
        this.nfe = this.nfe + 1;
        
        % Decode and Evaluate
        if 1
            [z, sol] = this.problem.func(xhat);
        else
            z = 0;
            sol = [];
        end
        end
        
% Iterations
function pso = iterate(pso)

% Decision Vector Size
var_size = pso.problem.dim;

% Move Particles
for i = 1:pso.pop_size
    
    % Update Velocity
    pso.pop(i).velocity = pso.w * pso.pop(i).velocity ...
        + pso.c1 * rand(1,var_size).*(pso.pop(i).best.position - pso.pop(i).position) ...
        + pso.c2 * rand(1,var_size).*(pso.best_sol.position - pso.pop(i).position);
    
    % Update Position
    pso.pop(i).position = pso.pop(i).position + pso.pop(i).velocity;
    
    % Mirror Velocities in Case of Limit Violation
    mirror_flag = (pso.pop(i).position < 0) | (pso.pop(i).position > 1);
    pso.pop(i).velocity(mirror_flag) = -pso.pop(i).velocity(mirror_flag);
    
    % Evaluate New Solution
    pso.pop(i) = eval(pso,pso.pop(i));
    
    % Comapre to Personal Best
    if pso.pop(i).obj_value < pso.pop(i).best.obj_value
       
        % Update Personal Best
        pso.pop(i).best.position = pso.pop(i).position;
        pso.pop(i).best.solution = pso.pop(i).solution;
        pso.pop(i).best.obj_value = pso.pop(i).obj_value;
        
        % Compare to Global Best
        t = pso.pop(i).best.obj_value<pso.best_sol.obj_value;
        if  (pso.pop(i).best.obj_value<pso.best_sol.obj_value)
            pso.best_sol = pso.pop(i).best;
        end
    end
    pso.best_obj_value = pso.best_sol;
    
end

% Damp Inertia Weight
pso.w = pso.wdamp * pso.w;

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
           
          opt.iter = 0;
          opt.nfe = 0;
          opt.nfe_history = nan(1, max_iter);
          opt.best_obj_value_history = nan(1, max_iter);
          opt.run_time = 0;
          opt.avg_eval_time = 0;
          opt.must_stop = false;
          opt.best_obj_value = inf
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
