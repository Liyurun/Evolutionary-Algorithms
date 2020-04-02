function position = covariance_matrix_adaptation_evolution_strategy(obj, problem)

    opt = struct('max_iter', 10,...
                'mu',50,...
                'lambda' , 200,...
                'disp',1,...
                'problem',problem);
        cmaes = struct('name','CMAES');
        cmaes = opt;
        cmaes = start(cmaes);
        
            if cmaes.disp
                disp(['CMAES  ' 'started ...']);
                disp('Initializing population.');
            end       
        
        cmaes = initialize(cmaes);
         
        % Iterations (Main Loop)
                for it = 1:cmaes.max_iter
                    
                    % Set Iteration Counter
                    cmaes.iter = it;
                    
                    % Iteration Started
                    cmaes = iterate(cmaes);
                    
                    % Update Histories
                    cmaes.nfe_history(it) = cmaes.nfe;
                    cmaes.best_obj_value_history(it) = cmaes.best_sol.obj_value;
                    position = cmaes.best_sol.position;
                    
                    % Display Iteration Information
                    if cmaes.disp
                        disp(['Iteration ' num2str(cmaes.iter) ...
                              ': Best this. Value = ' num2str(cmaes.best_sol.obj_value) ...
                              ', position = ' num2str(position)]);
                    end
                     
                    
                    
                end
                

        
        
        
end

function cmaes = start(cmaes)

cmaes.empty_individual.position = [];
cmaes.empty_individual.obj_value = [];
cmaes.empty_individual.solution = [];
cmaes.pop= [];
cmaes.best_sol= [];
cmaes.iter= 0;
cmaes.nfe= 0;
cmaes.run_time= 0;
cmaes.avg_eval_time= 0;
cmaes.nfe_history= [];
cmaes.best_obj_value_history= [];
cmaes.must_stop= 0;
% reset

            cmaes.iter = 0;
            cmaes.nfe = 0;
            cmaes.nfe_history = nan(1, cmaes.max_iter);
            cmaes.best_obj_value_history = nan(1, cmaes.max_iter);
            cmaes.run_time = 0;
            cmaes.avg_eval_time = 0;

end

% Iterations
function cmaes = iterate(cmaes)
% Iteration Counter
            g = cmaes.iter;
            
            % Decision Vector Size
            var_size = [1 cmaes.problem.dim];
            
            % Load Parameters
            w = cmaes.params.w;
            mu_eff = cmaes.params.mu_eff;
            cs = cmaes.params.cs;
            ds = cmaes.params.ds;
            ENN = cmaes.params.ENN;
            cc = cmaes.params.cc;
            c1 = cmaes.params.c1;
            cmu = cmaes.params.cmu;
            hth = cmaes.params.hth;
            M = cmaes.params.M;
            C = cmaes.params.C;
            pc = cmaes.params.pc;
            ps = cmaes.params.ps;
            sigma = cmaes.params.sigma;
            
            % Generate Samples (New Solutinos)
            cmaes.pop = repmat(cmaes.empty_individual, cmaes.lambda, 1);
            for i = 1:cmaes.lambda
                
                cmaes.pop(i).step = randn(var_size)*chol(C{g});
                cmaes.pop(i).position = M(g).position + sigma(g)*cmaes.pop(i).step;
                cmaes.pop(i).obj_value = cmaes.problem.func(cmaes.pop(i).position);
                
                if cmaes.pop(i).obj_value < cmaes.best_sol.obj_value
                    cmaes.best_sol = cmaes.pop(i);
                end
                
            end
            
            % Sort Population

            [obj_values, sort_order] = sort([cmaes.pop.obj_value]);
            cmaes.pop = cmaes.pop(sort_order);
            
            if g < cmaes.max_iter
                
                % Update Mean
                M(g+1).step = 0;
                for j = 1:cmaes.mu
                    M(g+1).step = M(g+1).step + w(j)*cmaes.pop(j).step;
                end
                M(g+1).position = M(g).position + sigma(g)*M(g+1).step;
                M(g+1).obj_value = cmaes.problem.func(M(g+1).position);
                if  M(g+1).obj_value < cmaes.best_sol.obj_value
                    cmaes.best_sol = M(g+1);
                end
                
                % Update Step Size
                ps{g+1} = (1-cs)*ps{g} + sqrt(cs*(2-cs)*mu_eff)*M(g+1).step/chol(C{g})';
                sigma(g+1) = sigma(g)*exp(cs/ds*(norm(ps{g+1})/ENN-1))^0.3;
                
                % Update Covariance Matrix
                if norm(ps{g+1})/sqrt(1-(1-cs)^(2*(g+1))) < hth
                    hs=1;
                else
                    hs=0;
                end
                delta = (1-hs)*cc*(2-cc);
                pc{g+1} = (1-cc)*pc{g} + hs*sqrt(cc*(2-cc)*mu_eff)*M(g+1).step;
                C{g+1} = (1-c1-cmu)*C{g} + c1*(pc{g+1}'*pc{g+1} + delta*C{g});
                for j = 1:cmaes.mu
                    C{g+1} = C{g+1} + cmu*w(j)*cmaes.pop(j).step'*cmaes.pop(j).step;
                end
                
                % If Covariance Matrix is not Positive Defenite or Near Singular
                [V, E] = eig(C{g+1});
                if any(diag(E)<0)
                    E = max(E,0);
                    C{g+1}=V*E/V;
                end
                
            end
            
            % Store Parameters
            cmaes.params.M = M;
            cmaes.params.C = C;
            cmaes.params.pc = pc;
            cmaes.params.ps = ps;
            cmaes.params.sigma = sigma;
end

function cmaes = initialize(cmaes)

            % Decision Variables Count
            var_count = cmaes.problem.dim;
            
            % Decision Vector Size
            var_size = [1 var_count];
            
            % Parent Weights
            w = log(cmaes.mu + 0.5) - log(1:cmaes.mu);
            w = w/sum(w);
            
            % Number of Effective Solutions
            mu_eff = 1/sum(w.^2);
            
            % Step Size Control Parameters (c_sigma and d_sigma);
            sigma0 = 0.3;
            cs = (mu_eff + 2)/(var_count + mu_eff + 5);
            ds = 1 + cs + 2*max(sqrt((mu_eff-1)/(var_count+1))-1,0);
            ENN = sqrt(var_count)*(1-1/(4*var_count)+1/(21*var_count^2));
            
            % Covariance Update Parameters
            cc = (4 + mu_eff/var_count)/(4 + var_count + 2*mu_eff/var_count);
            c1 = 2/((var_count+1.3)^2+mu_eff);
            alpha_mu = 2;
            cmu = min(1-c1,alpha_mu*(mu_eff-2+1/mu_eff)/((var_count+2)^2+alpha_mu*mu_eff/2));
            hth = (1.4+2/(var_count+1))*ENN;
            
            % Store Params
            cmaes.params.w = w;
            cmaes.params.mu_eff = mu_eff;
            cmaes.params.sigma0 = sigma0;
            cmaes.params.cs = cs;
            cmaes.params.ds = ds;
            cmaes.params.ENN = ENN;
            cmaes.params.cc = cc;
            cmaes.params.c1 = c1;
            cmaes.params.alpha_mu = alpha_mu;
            cmaes.params.cmu = cmu;
            cmaes.params.hth = hth;
            
            % Initialize Step Sizes for Individuals
            cmaes.empty_individual.step = zeros(var_size);
            
            % Initialize first Mean (M)
            cmaes.params.M = repmat(cmaes.empty_individual, cmaes.max_iter, 1);
            x = rand([1 cmaes.problem.dim]);
            cmaes.params.M(1) = new_individual(cmaes,x);
            
            % Initialize first Covariance Matrix (C)
            cmaes.params.C = cell(cmaes.max_iter, 1);
            cmaes.params.C{1} = eye(var_count);
            
            % Initialize first pc
            cmaes.params.pc = cell(cmaes.max_iter, 1);
            cmaes.params.pc{1} = zeros(var_size);
            
            % Initialize first ps
            cmaes.params.ps = cell(cmaes.max_iter, 1);
            cmaes.params.ps{1} = zeros(var_size);
            
            % Initialize first Step Size (sigma)
            cmaes.params.sigma = zeros(cmaes.max_iter, 1);
            cmaes.params.sigma(1) = sigma0;
            
            % Initialize Population
            cmaes.pop = [];
            
            % Initialize Best Solution Ever Found
            cmaes.best_sol = cmaes.params.M(1);            
            
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


        