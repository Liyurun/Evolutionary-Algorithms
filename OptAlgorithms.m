
classdef OptAlgorithms < handle
% =============================================================================
% This is the class of the functions of optimization algorithms.
%
% =============================================================================   

    properties (Constant)
        FUNC        = @func;   % the objective function
        swarm       = 2;               % the population amount
        sample      = 10;             % the loop count for evolution
        dataPoint   = 20;              % the amount of data observation
        DE_strategy = 5;                % the strategy of DE kernel
    end


% Lower level algorithms, in charge of continuous decision variables optimization 
% -----------------------------------------------------------------------------
%   DE
    methods (Static = true, Access = 'public')

        function position = Differential_Evolution(obj, problem)
            position = differential_evolution(obj, problem);

    end % DE
    end

%   PSO    
    methods (Static = true, Access = 'public')

        function position = Particle_Swarm_Optimization(obj, problem)
% -----------------------------------------------------------------------------
% Particle Swarm Optimization algorithm (PSO)
%
%  PSO optimizes a problem by having a population of candidates
%  (particles), and moving these particles around in the search space
%  according to mathematical formula ovet the particle's position and
%  velocity. Each particle's movement is influenced by its own local
%  best-known position but, is also guided toward the best-known positions
%  in the search space, which are updated as better positions are found by
%  other particles
%
%  Parameter:
%        - params. It is the specified parameters from the main function.
%        And let the main function informed that which parameters are need
%        to be optimized
% -----------------------------------------------------------------------------

            position = particle_swarm_optimization(obj, problem);
            
        end
    end


%   GA    
    methods (Static = true, Access = 'public')

        function position = Genetic_Algorithm(obj, problem)

            position = genetic_algorithm(obj, problem);
            
        end
    end
    
    %   CMAES   
    methods (Static = true, Access = 'public')

        function position = Covariance_Matrix_Adaptation_Evolution_Strategy(obj, problem)

            position =covariance_matrix_adaptation_evolution_strategy(obj, problem);
            
        end
    end
    % SA 
    methods (Static = true, Access = 'public')

        function position = Simulated_Annealing(obj, problem)

            position =simulated_annealing(obj, problem);
            
        end
    end
     % ABC
    methods (Static = true, Access = 'public')

        function position = Artificial_Bee_Colony(obj, problem)

            position =artificial_bee_colony(obj, problem);
            
        end
    end
         % ACO
    methods (Static = true, Access = 'public')

        function position =    Ant_Colony_Optimization(obj, problem)

            position = ant_colony_optimization(obj, problem);
            
        end
    end
             % BA
    methods (Static = true, Access = 'public')

        function position =    Bees_Algorithm(obj, problem)

            position = bees_algorithm(obj, problem);
            
        end
    end
      % BBO
    methods (Static = true, Access = 'public')

        function position = Biogeography_Based_Optimization(obj, problem)

            position = biogeography_based_optimization(obj, problem);
            
        end
    end
      % FA
    methods (Static = true, Access = 'public')

        function position = Firefly_Algorithm(obj, problem)

            position = firefly_algorithm(obj, problem);
            
        end
    end

          % HS
    methods (Static = true, Access = 'public')

        function position = Harmony_Search(obj, problem)

            position = harmony_search(obj, problem);
            
        end
    end
    
              % ICA
    methods (Static = true, Access = 'public')

        function position = Imperialist_Competitive_Algorithm(obj, problem)

            position = imperialist_competitive_algorithm(obj, problem);
            
        end
    end
    
    
   % IWO
    methods (Static = true, Access = 'public')

        function position = Invasive_Weed_Optimization(obj, problem)

            position = invasive_weed_optimization(obj, problem);
            
        end
    end
    
    
    % TLBO
    methods (Static = true, Access = 'public')

        function position = Teaching_Learning_Based_Optimization(obj, problem)

            position = teaching_learning_based_optimization(obj, problem);
            
        end
    end
    
    
  
    
    
    
    
%	MCMC
    methods (Static = true, Access = 'public')

        function [xValue, yValue] = Markov_Chain_Monte_Carlo(obj, params)
%------------------------------------------------------------------------------
% Markov chain Monte Carlo (MCMC)simulation
%
% The central problem is that of determining the posterior probability for
% parameters given the data, p(thelta|y). With uninformative prior distribution
% and normalization constant p(y), the task reduced to that of maximizing the
% likelihood of data, p(y|thelta).
%
% "Delayed rejection" was implemented in order to increase the acceptance ratio.
%
% If the Jacobian matrix can be obtained, it will be used to generate R
%------------------------------------------------------------------------------


            startTime = clock;

            if nargin < 2
                error('OptAlgorithms.MCMC: There are no enough inputs \n');
            end

            % Get the MCMC options
            opt = OptAlgorithms.getOptions_MCMC(obj, params);

            % Inilization of the guessing point
            [oldpar, SS, Jac] = OptAlgorithms.initPoint(opt);

            % Preallocation
            lasti = 0; chaincov = []; chainmean = []; wsum = [];
            accepted =0;
            n0 = 1;

            chain = zeros(opt.nsamples, opt.Nparams+1);

            % Construct the error standard deviation: sigma square
            sumSquare = SS;
            if SS < 0
                sumSquare = exp(SS);
            end
            sigmaSqu = sumSquare / (opt.nDataPoint - opt.Nparams); 
            sigmaSqu_0 = sigmaSqu;

            if ~isempty(Jac)

                % Construct the R = chol(covMatrix) for candidate generatio
                [~, S, V] = svd(Jac);

                % Truncated SVD
                S = diag(S);
                S(S < 1e-3) = 1e-3;

                % Fisher information matrix
                % covariance matrix = sigmaSqu * inv(Fisher information matrix) = v'*s^{-2}*v
                covMatrix = V * diag(1./S.^2)* V' * sigmaSqu;

                % Construct the R = chol(covMatrix) for candidate generation
                R = chol(covMatrix .* 2.4^2 ./ opt.Nparams);

            else

                % If the Jacobian matrix cannot be obtained,
                % a set of samples is used to generate R
                R = OptAlgorithms.burnInSamples(opt);

            end

%------------------------------------------------------------------------------

% Main loop
%------------------------------------------------------------------------------
            for j = 1:(opt.nsamples + opt.burn_in)

                if j > opt.burn_in
                    fprintf('Iter: %4d -------- Accept_ratio: %3d%% ---------- Minimum: %g ---------- \n',...
                        j, fix(accepted / j * 100), SS);
                    fprintf('%10.3g | ', oldpar); fprintf('\n');
                end

                accept = false;

                % Generate the new proposal point with R
                newpar = oldpar + randn(1,opt.Nparams) * R;

                % Check the boundary limiation
                col = find(newpar < opt.bounds(1, :));
                newpar(col) = opt.bounds(1, col);
                col = find(newpar > opt.bounds(2, :));
                newpar(col) = opt.bounds(1, col);

                % Calculate the objective value of the new proposal
                newSS = feval( OptAlgorithms.FUNC, newpar );

                % The Metropolis probability
                rho12 = exp( -0.5 * (newSS - SS) / sigmaSqu);

                % The new proposal is accepted with Metropolis probability
                if rand <= min(1, rho12)
                    accept   = true;
                    oldpar   = newpar;
                    SS       = newSS;
                    accepted = accepted + 1;
                end

                % If the poposal is denied, a Delayed Rejection procedure is adopted
                % in order to increase the acceptance ratio
                if ~accept && opt.DelayedRejection

                    % Shrink the searching domain by the factor 1/10
                    newpar2 = oldpar + randn(1, opt.Nparams) * (R ./ 10);

                    % Check the boundary limitation of the new generated point
                    col = find(newpar2 < opt.bounds(1, :));
                    newpar2(col) = opt.bounds(1, col);
                    col = find(newpar2 > opt.bounds(2, :));
                    newpar2(col) = opt.bounds(1, col);

                    % Calculate the objective value of the new proposal
                    newSS2 = feval( OptAlgorithms.FUNC, newpar2 );

                    % The conventional version of calculation
%                   rho32 = exp( -0.5 * (newSS - newSS2) / sigmaSqu);

%                   q2 = exp( -0.5 * (newSS2 - SS) / sigmaSqu);
%                   q1 = exp( -0.5 * (norm((newpar2 - newpar) * inv(R))^2 - norm((oldpar - newpar) * inv(R))^2));

%                   if rho32 == Inf
%                       rho13 = 0;
%                   else
%                       rho13 = q1 * q2 * (1-rho32) / (1-rho12);
%                   end

                    % The speed-up version of above calculation
                    q1q2 = exp( -0.5 * ((newSS2 - SS) / sigmaSqu + ...
                        ((newpar2 - newpar) * (R \ (R' \ (newpar2' - newpar')))...
                        - (oldpar - newpar) * (R \ (R' \ (oldpar' - newpar'))))));

                    rho13 = q1q2 * expm1(-0.5 * (newSS - newSS2) / sigmaSqu) / expm1(-0.5 * (newSS - SS) / sigmaSqu);

                    if rand <= min(1, rho13)
                        oldpar   = newpar2;
                        SS       = newSS2;
                        accepted = accepted + 1;
                    end

                end

                % During the burn-in period, if the acceptance rate is extremly high or low,
                % the R matrix is manually adjusted
                if j <= opt.burn_in
                    if mod(j, 50) == 0
                        if accepted/j < 0.05
                            fprintf('Acceptance ratio %3.2f smaller than 5 %%, scaled \n', accepted/j*100);
                            R = R ./ 5;
                        elseif accepted/j > 0.95
                            fprintf('Acceptance ratio %3.2f largeer than 95 %%, scaled \n', accepted/j*100);
                            R = R .* 5;
                        end
                    end
                end

                % After the burn-in period, the chain is stored
                if j > opt.burn_in
                    chain(j - opt.burn_in, 1:opt.Nparams)   = oldpar;
                    chain(j - opt.burn_in, opt.Nparams + 1) = SS;

                    temp = chain(j-opt.burn_in, :);
                    save('chainData.dat', 'temp', '-ascii', '-append');
                end


                % Updata the R according to previous chain
                % Check the convergence condition
                if mod((j-opt.burn_in), opt.convergInt) == 0 && (j-opt.burn_in) > 0

                    [chaincov, chainmean, wsum] = OptAlgorithms.covUpdate(chain((lasti+1):j-opt.burn_in, 1:opt.Nparams),...
                        1, chaincov, chainmean, wsum);

                    lasti =j;
                    R = chol(chaincov + eye(opt.Nparams)*1e-7);

                    criterion = OptAlgorithms.Geweke( chain(1:j-opt.burn_in, 1:opt.Nparams) );

                    if all( abs(criterion) < opt.criterionTol ) || j == opt.nsamples+opt.burn_in
                        maxIter = j;
                        break
                    end

                end


                % Updata the sigma^2 according to the current objective value
                if SS < 0, sumSquare = exp(SS); end
                sigmaSqu  = 1 / OptAlgorithms.GammarDistribution( 1, 1, (n0 + opt.nDataPoint)/2,...
                    2 / (n0 * sigmaSqu_0 + sumSquare));


            end
%------------------------------------------------------------------------------

% Post-process
%------------------------------------------------------------------------------

            clear chain;

            % Generate the population for figure plot
            Population = OptAlgorithms.convertASCIItoMAT_MCMC(maxIter, opt);

            OptAlgorithms.FigurePlot(Population, opt);

            [yValue, row] = min(Population(:, opt.Nparams+1));
            xValue        = Population(row, 1:opt.Nparams);

            % Gather some useful information and store them
            result.optTime         = etime(clock, startTime) / 3600;
            result.covMatrix       = R' * R;
            result.Iteration       = maxIter;
            result.criterion       = criterion;
            result.accepted        = fix(accepted/maxIter * 100);
            result.xValue          = xValue;
            result.yValue          = yValue;

            save(sprintf('result_%3d.mat', fix(rand*1000)), 'result');
            fprintf('Markov chain simulation finished and the result is saved as result.mat \n');

        end

        function opt = getOptions_MCMC(obj, params)
%------------------------------------------------------------------------------
% The parameters for the Optimizer
%
% Parameter:
%       - params. It is the specified parameters from the main function.
%       And let the main function informed that which parameters are needed
%       to be optimized.
%
% Return:
%       - opt.
%           + opt.Nparams. The number of optimized parameters
%           + opt.bounds. 2*nCol matrix of the parameter limitation
%           + opt.nsamples. The pre-defined maximal iteration
%           + opt.criterionTol. The error tolerance to stop the algorithm
%           + opt.burn_in. The burn-in period before adaptation begin
%           + opt.burnInSwarm. This is specific used in the generation of R matrix
%               when Jacobian information is not available. A sample of swarm is initialized
%           + opt.convergInt. The integer for checking of the convergence
%           + opt.rejectValue. This is specific used in the generation of R matrix
%               when Jacobian information is not available. In the generated sample,
%               the proposals whose objective value is larger than this will be rejected
%           + opt.nDataPoint. The number of observations in the measured data
%           + opt.DeylayedRejection. By default DR= 1
%           + opt.Jacobian. If Jacobian matrix is available, set it to true
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgorithms.getOptions_MCMC: There are no enough input arguments \n');
            end

            opt = [];

            opt.Nparams       = length(fieldnames(params));
            opt.bounds        = obj.paramBound';

            opt.nsamples      = OptAlgorithms.sample;
            opt.criterionTol  = 0.0001;
            opt.burn_in       = 100;
            opt.burnInSwarm   = 5000;
            opt.convergInt    = 50;
            opt.rejectValue   = 10000;
            opt.nDataPoint    = OptAlgorithms.dataPoint;

            opt.Jacobian      = false;
            opt.DelayedRejection = true;

            [row, col] = size(opt.Nparams);
            if row > 1 || col > 1
                error('OptAlgorithms.getOptions_MCMC: The initialized dimension of the parameter set might be wrong \n');
            end

            [row, col] = size(opt.bounds);
            if col ~= opt.Nparams || row ~= 2
                error('OptAlgorithms.getOptions_MCMC: Please check your setup of the range of parameters \n');
            end

        end

        function [oldpar, SS, Jac]= initPoint(opt)
%------------------------------------------------------------------------------
% Generate the initially guessing point for the MCMC algorithm
%------------------------------------------------------------------------------


            oldpar = rand(1,opt.Nparams);

            oldpar(:, 1:opt.Nparams) = opt.bounds(1,:) + ...
                oldpar(:, 1:opt.Nparams) .* (opt.bounds(2,:) - opt.bounds(1,:));

            % Check the boundary limiation
%            col = find(oldpar < opt.bounds(1, :));
%            oldpar(col) = opt.bounds(1, col);
%            col = find(oldpar > opt.bounds(2, :));
%            oldpar(col) = opt.bounds(1, col);

            % Get the resudual value and Jacobian matrix of the guessing point
            if opt.Jacobian
                [SS, Jac] = feval( OptAlgorithms.FUNC, oldpar );
            else
                SS = feval( OptAlgorithms.FUNC, oldpar );
                Jac = [];
            end

        end

        function R = burnInSamples(opt)
%------------------------------------------------------------------------------
% The routine that is used for generating samples in cases that Jocabian matrix
% is not available
%------------------------------------------------------------------------------


            Swarm = OptAlgorithms.swarm;

            % Generate random sample with swarm scale
            ParSwarm = rand(Swarm, opt.Nparams+1);

            % Keep all the swarm in the searching domain
            ParSwarm(:, 1:opt.Nparams) = repmat(opt.bounds(1,:), Swarm,1) + ...
                ParSwarm(:, 1:opt.Nparams) .* repmat( (opt.bounds(2,:) - opt.bounds(1,:)), Swarm,1 );

            % All the proposals are simulated
            ParSwarm(:, opt.Nparams+1) = arrayfun( @(idx) feval(OptAlgorithms.FUNC,...
                ParSwarm(idx, 1:opt.Nparams) ), 1:Swarm);

            % If the objective value is too big in the case, the proposals are deleted
            row = find(ParSwarm(:, opt.Nparams+1) > opt.rejectValue);
            ParSwarm(row, :) = [];

            % Calculate the covariance matrix
            [chaincov, ~, ~] = OptAlgorithms.covUpdate(ParSwarm(:, 1:opt.Nparams), 1, [], [], []);

            % Cholesky decomposition
            R = chol( chaincov + eye(opt.Nparams) * 1e-7 );

        end

        function z = Geweke(chain, a, b)
%------------------------------------------------------------------------------
% Geweke's MCMC convergence diagnostic
%------------------------------------------------------------------------------


            if nargin < 3
                a = 0.5;
                if nargin < 2
                    b = 0.8;
                end
            end

            [n, par] = size(chain);
            na = floor(a*n);
            nb = floor(b*n);

            if (n-na) / n >= 1
                error('OptAlgorithms.Geweke: Error with na and nb \n');
            end

            m1 = mean(chain(na:(nb-1),:));
            m2 = mean(chain(nb:end,:));

            % Spectral estimates for variance
            sa = OptAlgorithms.spectrum0(chain(na:(nb-1),:));
            sb = OptAlgorithms.spectrum0(chain(nb:end,:));

            z = (m1 - m2) ./ (sqrt( sa / (nb-na) + sb / (n-nb+1)));

        end

        function s = spectrum0(x)
%------------------------------------------------------------------------------
% Spectral density at frequency zero
%------------------------------------------------------------------------------


            [m, n] = size(x);
            s = zeros(1,n);

            for i = 1:n
                spec = OptAlgorithms.spectrum(x(:,i),m);
                s(i) = spec(1);
            end

        end

        function [y, f] = spectrum(x, nfft, nw)
%------------------------------------------------------------------------------
% Power spectral density using Hanning window
%------------------------------------------------------------------------------


            if nargin < 3 || isempty(nw)
                nw = fix(nfft/4);
                if nargin < 2 || isempty(nfft)
                    nfft = min(length(x),256);
                end
            end

            noverlap = fix(nw/2);

            % Hanning window
            w = .5*(1 - cos(2*pi*(1:nw)' / (nw+1)));
            % Daniel
%            w = [0.5;ones(nw-2,1);0.5];
            n = length(x);

            if n < nw
                x(nw) = 0;
                n = nw;
            end

            x = x(:);

            k = fix((n - noverlap) / (nw - noverlap)); % no of windows
            index = 1:nw;
            kmu = k * norm(w)^2; % Normalizing scale factor
            y = zeros(nfft,1);

            for i = 1:k
%                xw = w.*detrend(x(index),'linear');
                xw = w .* x(index);
                index = index + (nw - noverlap);
                Xx = abs(fft(xw,nfft)).^2;
                y = y + Xx;
            end

            y  = y * (1 / kmu); % normalize

            n2 = floor(nfft / 2);
            y  = y(1:n2);
            f  = 1 ./ n * (0:(n2-1));

        end

        function [xcov, xmean, wsum] = covUpdate(x, w, oldcov, oldmean, oldwsum)
%------------------------------------------------------------------------------
% Recursive update the covariance matrix
%------------------------------------------------------------------------------


            [n, p] = size(x);

            if n == 0
                xcov = oldcov;
                xmean = oldmean;
                wsum = oldwsum;
                return
            end

            if nargin < 2 || isempty(w)
                w = 1;
            end

            if length(w) == 1
                w = ones(n,1) * w;
            end

            if nargin > 2 && ~isempty(oldcov)

                for i = 1:n
                    xi     = x(i,:);
                    wsum   = w(i);
                    xmeann = xi;
                    xmean  = oldmean + wsum / (wsum + oldwsum) * (xmeann - oldmean);

                    xcov =  oldcov + wsum ./ (wsum + oldwsum - 1) .* (oldwsum / (wsum + oldwsum) ...
                        .* ((xi - oldmean)' * (xi - oldmean)) - oldcov);
                    wsum    = wsum + oldwsum;
                    oldcov  = xcov;
                    oldmean = xmean;
                    oldwsum = wsum;
                end

            else

                wsum  = sum(w);
                xmean = zeros(1,p);
                xcov  = zeros(p,p);

                for i = 1:p
                    xmean(i) = sum(x(:,i) .* w) ./ wsum;
                end

                if wsum > 1
                    for i = 1:p
                        for j = 1:i
                            xcov(i,j) = (x(:,i) - xmean(i))' * ((x(:,j) - xmean(j)) .* w) ./ (wsum - 1);
                            if (i ~= j)
                                xcov(j,i) = xcov(i,j);
                            end
                        end
                    end
                end

            end

        end

        function Population = convertASCIItoMAT_MCMC(maxIter, opt)
%------------------------------------------------------------------------------
% Load and convert the data in ascii format to the mat format
% then generate the population for statistics
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgorithms.convertASCIItoMAT: There are no enough input arguments \n');
            end

            load('chainData.dat');

            % Discard the former 50% chain, retain only the last 50%
            if maxIter < opt.nsamples + opt.burn_in
                idx = floor(0.5 * maxIter);
            else
                idx = floor(0.5 * (opt.nsamples + opt.burn_in));
            end

            eval(sprintf('chainData(1:idx, :) = [];'));
            Population = chainData;
            
            save('population.dat', 'Population', '-ascii');

        end

	end % MCMC


%   PRML
    methods (Static = true, Access = 'public')

        function [xValue, yValue] = Parallel_Riemann_Metropolis_Adjusted_Langevin(obj, params)
%------------------------------------------------------------------------------ 
% Riemannian manifold Metropolis adjusted Langevin with parallel tempering (PRML)
%
% The Langevin algorithm is combined withe Metropolis probability, which leads to
% the MALA method in the literatures.
% Riemann manifold is introduced into the MALA method to enhance the searching capacity
% However, it still has the possibility to trap into the local optima, so the parallel
% tempering technique is employed to tackle the multi-variable optmization
%
% Parameter:
%       - obj. Some options for the algorithms that is from the main function
%       - params. It is the specified parameters from the main function.
%
% Returns:
%       - xValue. The vector of the optimal parameter
%       - yValue. The value of the objective function
%------------------------------------------------------------------------------


            global accepted;

            startTime = clock;

            if nargin < 2
                error('OptAlgorithms.PRML: There are no enough input arguments \n');
            end

            % Get the PRML options
            opt = OptAlgorithms.getOptions_PRML(obj, params);

            % Preallocation
            accepted = 0; n0 = 1;

            % The temperature matrix of the parallel tempering
            Temperatures = zeros(opt.nsamples, opt.Nchain);
            Temperatures(:, 1) = ones(opt.nsamples, 1);
            Temperatures(1, :) = opt.temperature;

            % Initialization 
            [states, MetricTensor] = OptAlgorithms.initChainP((1./opt.temperature), opt);

            chain = zeros(opt.Nchain, opt.Nparams+1, opt.nsamples);
%            sigmaChain = zeros(opt.nsamples, opt.Nchain);

            % Calculate the initial sigma square values
            sumSquare = states(:, opt.Nparams+1);
            if any(sumSquare < 0)
                row = find(sumSquare < 0);
                sumSquare(row) = exp( sumSquare(row) );
            end
            sigmaSqu  = sumSquare/ (opt.nDataPoint - opt.Nparams);

            sigmaSqu0 = sigmaSqu;

%------------------------------------------------------------------------------

% Main loop
%------------------------------------------------------------------------------
            for i = 1:opt.nsamples

                if i > opt.burn_in && mod(i, opt.convergInt) == 0
                    fprintf('Iter: %3d ----- Accept_ratio: %2d%% ', i, fix(accepted/i*100));

                    % Abstract best information so far from the population and display it
                    [minValue, row] = min(states(1:opt.Nchain, opt.Nparams+1));
                    fprintf('----------------  Minimum: %3g  ---------------- \n', minValue);
                    fprintf('%10.3g | ', states(row, 1:opt.Nparams) ); fprintf('\n');
                end

                % PRML: Evolution of the chains
                [states, MetricTensor] = OptAlgorithms.PRML_sampler(states, MetricTensor, sigmaSqu, 1./Temperatures(i,:), opt);

                % Store the chain after burn-in period
                if i > opt.burn_in
                    chain(:,:,i) = states;

                    for j = 1:opt.Nchain
                        temp = states(j, :);
                        save(sprintf('chainData_%d.dat', j), 'temp', '-ascii', '-append');
                    end
                end

                % Implement the chain swap
                if mod(i, opt.swapInt) == 0
                    states = OptAlgorithms.chainSwap(states, sigmaSqu, 1./Temperatures(i,:), opt);
                end

                % In each opt.convergInt interval, check the convergence diagnostics
                if mod(i, opt.convergInt) == 0 && i > opt.convergInt

                    criterion = OptAlgorithms.GelmanR_statistic(i, chain(:,:,1:i), opt);

                    if all(criterion < opt.criterionTol) || i == opt.nsamples
                        maxIter = i;
                        break
                    end

                end

                % Variance of error distribution (sigma) was treated as a parameter to be estimated.
                sumSquare = states(:, opt.Nparams+1);
                if any(sumSquare < 0)
                    row = find(sumSquare < 0);
                    sumSquare(row) = exp( sumSquare(row) );
                end

                for k = 1:opt.Nchain
                    sigmaSqu(k)  = 1 ./ OptAlgorithms.GammarDistribution(1, 1, (n0 + opt.nDataPoint)/2, ...
                        2 / (n0 * sigmaSqu0(k) + sumSquare(k)));
%                    sigmaChain(i,k) = sigmaSqu(k)';
                end

                % Temperature dynamics
                Temperatures = OptAlgorithms.temperatureDynamics(i, states, Temperatures, sigmaSqu, opt);


            end  % for i = 1:opt.nsamples
%------------------------------------------------------------------------------

% Post-process
%------------------------------------------------------------------------------

            clear chain;

            % Generate the population for figure plot
            Population = OptAlgorithms.convertASCIItoMAT_PRML(maxIter, opt);

            % Plot the population
            OptAlgorithms.FigurePlot(Population, opt)

            [yValue, row] = min(Population(:,opt.Nparams+1));
            xValue        = Population(row,1:opt.Nparams);

            % Gather some useful information and store them
            result.optTime        = etime(clock, startTime)/3600;
            result.convergDiago   = criterion;
            result.NChain         = opt.Nchain;
            result.Ieteration     = maxIter;
            result.Temperatures   = Temperatures(1:opt.nsamples, :);
            result.acceptanceRate = fix(accepted/maxIter*100);
            result.correlation    = corrcoef(Population(1:opt.Nparams));
            result.population     = Population;
            result.xValue         = xValue;
            result.yValue         = yValue;

            save(sprintf('result_%3d.mat', fix(rand*1000)), 'result');
            fprintf('Markov chain simulation finished and the result saved as result.mat \n');

        end

        function opt = getOptions_PRML(obj, params)
%------------------------------------------------------------------------------
% The parameters for the optimizer
%
% Parameter:
%       - obj. The parameter options from the main function
%       - params. It is the specified parameters from the main function
% 
% Return:
%       - opt.
%           + temperature. The vector of the temperature of the parallel tempering
%           + Nchain. The number of the candidates
%           + Nparams. The number of the optimized parameters
%           + nsamples. The length of the sampling chain
%           + bounds. The boundary limitation of the sampling
%           + burn_in. The burn-in period that is discarded
%           + convergInt. The interval to check the convergence criterion 
%           + swapInt. The interval to swap the N chains when sampling
%           + nDataPoint. The data point of the y
%           + criterionTol. The stopping tolerance
%           + epsilon. The epsilon value in the PRML proposal formula
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgorithms.getOptions_PRML: There are no enough input arguments \n');
            end

            opt = [];

            opt.temperature       = [1 10 100 500];
            opt.Nchain            = length(opt.temperature);
            opt.Nparams           = length(fieldnames(params));
            opt.nsamples          = OptAlgorithms.sample;
            opt.bounds            = obj.paramBound';

            opt.burn_in           = 100;
            opt.convergInt        = 50;
            opt.swapInt           = 100;
            opt.nDataPoint        = OptAlgorithms.dataPoint;
            opt.criterionTol      = 1.1;
            opt.epsilon           = 0.15;

            if mod(opt.nsamples, opt.convergInt) ~= 0
                error('OptAlgorithms.getOptions_PRML: Please set the samples be devisible to the convergInt \n');
            end

            [row, col] = size(opt.Nparams);
            if row > 1 || col > 1
                error('OptAlgorithms.getOptions_PRML: The initialized dimension of the parameter set might be wrong \n');
            end

            [row, col] = size(opt.bounds);
            if col ~= opt.Nparams || row ~= 2
                error('OptAlgorithms.getOptions_PRML: Please check your setup of the range of parameters \n');
            end

        end

        function [initChain, MetricTensor] = initChainP(Beta, opt)
%------------------------------------------------------------------------------
% Generate initial chains for the PRML algorithm
%
% parameter:
%       - Beta. The inverse of the temperature vector
%       - opt. The options of the PRML
%
% Return:
%       - initChain. The Nchain x (Nparams+1) matrix
%       - MetricTensor.
%           + G. The Fisher information matrix (Metric tensor under Riemann manifold)
%           + GradL. The gradient of the log-density distribution L.
%           + sqrtInvG. The inverse of the Fisher information matrix
%------------------------------------------------------------------------------


            % Initilization of the chains
            initChain = rand(opt.Nchain, opt.Nparams+1);
            MetricTensor = cell(1, opt.Nchain);

            % Use vectorization to speed up. Keep the parameters in the domain
            initChain(:, 1:opt.Nparams) = repmat(opt.bounds(1, :), opt.Nchain, 1) + ...
                initChain(:, 1:opt.Nparams) .* repmat( (opt.bounds(2,:)-opt.bounds(1,:)), opt.Nchain, 1 );

            for j = 1:opt.Nchain

                % Simulation
                [Res, Jac] = feval( OptAlgorithms.FUNC, initChain(j, 1:opt.Nparams) );

                initChain(j, opt.Nparams+1) = Res' * Res;
                SigmaSqu = (initChain(j, opt.Nparams+1) / (opt.nDataPoint - opt.Nparams));

                % Prepare the information for the PRML proposal kernel
                %   metric tensor G, gradient vector, square root of inverse G
                MetricTensor{j}.G     = Beta(j) .* (Jac' * (1/SigmaSqu) * Jac);
                MetricTensor{j}.GradL = - Jac' * Res / SigmaSqu;
                if rank(MetricTensor{j}.G) ~= opt.Nparams
                    invG = pinv(MetricTensor{j}.G + eye(opt.Nparams)*1e-10);
                else
                    invG = inv(MetricTensor{j}.G);
                end

                % If inv(G) is not sysmmetric, SVD is used to approximate
                [R, p] = chol(invG);
                if p == 0
                    MetricTensor{j}.sqrtInvG = R;
                else
                    [U,S,V] = svd( invG );
                    S = diag(S); S(S<1e-10) = 1e-10;
                    MetricTensor{j}.sqrtInvG = U * diag(sqrt(S)) * V';
                end

            end

        end

        function [states, MetricTensor] = PRML_sampler(states, MetricTensor, sigmaSqu, Beta, opt)
%------------------------------------------------------------------------------
% The Metropolis adjusted Langevin algorithms is used to generate new proposal
%
% Parameter:
%       - states. The Nchain x (Nparams+1) matrix
%       - MetricTensor. A struct that contains the Fisher information, gradient
%           of the log-density distribution, inverse of the metric tensor
%       - sigmaSqu. The sigma square matrix (the covariance matrix)
%       - Beta. Inverse of the temperature vector in parallel tempering
%       - opt. The PRML options
%
% Return:
%       - states. Same
%       - MetricTensor. Same
%------------------------------------------------------------------------------


            global accepted;

            for j = 1: opt.Nchain

                % New proposal formula based on the Metropolis adjusted Langevin algorithm
                proposal = states(j,1:opt.Nparams) + 0.5 * opt.epsilon^2 * (MetricTensor{j}.G \...
                    MetricTensor{j}.GradL)'+ opt.epsilon * randn(1,opt.Nparams) * MetricTensor{j}.sqrtInvG;

                % Check the boundary limiation
                col = find(proposal < opt.bounds(1, :));
                proposal(col) = opt.bounds(1, col);
                col = find(proposal > opt.bounds(2, :));
                proposal(col) = opt.bounds(1, col);

                % Simulation of the new proposal
                [newRes, newJac] = feval( OptAlgorithms.FUNC, proposal );

                newSS = newRes' * newRes;
                SS    = states(j,opt.Nparams+1);

                % The Metropolis probability
                rho = (exp( -0.5*(newSS - SS) / sigmaSqu(j)))^Beta(j);

                % If the proposal is accepted
                if rand <= min(1, rho)

                    states(j, 1:opt.Nparams) = proposal;
                    states(j, opt.Nparams+1) = newSS;

                    % Prepare the information for the PRML proposal kernel
                    %   metric tensor G, gradient vector, square root of inverse G
                    MetricTensor{j}.G     = Beta(j) .* (newJac' * (1/sigmaSqu(j)) * newJac);
                    MetricTensor{j}.GradL = -newJac' * newRes / sigmaSqu(j);
                    if rank(MetricTensor{j}.G) ~= opt.Nparams
                        invG = pinv(MetricTensor{j}.G + eye(opt.Nparams)*1e-10);
                    else
                        invG = inv(MetricTensor{j}.G);
                    end

                    % If inv(G) is not sysmmetric, SVD is used to approximate
                    [R, p] = chol(invG);
                    if p == 0
                        MetricTensor{j}.sqrtInvG = R;
                    else
                        [U,S,V] = svd( invG );
                        S = diag(S); S(S<1e-10) = 1e-10;
                        MetricTensor{j}.sqrtInvG = U * diag(sqrt(S)) * V';
                    end

                    if j == 1, accepted = accepted + 1; end

                end

            end


        end

        function states = chainSwap(states, sigmaSqu, Beta, opt)
%------------------------------------------------------------------------------
% Swap of the chains at pre-determined intervals
%
% Parameter:
%       - states. The Nchain x (Nparams+1) matrix
%       - sigmaSqu. The sigma square matrix (the covariance matrix)
%       - Beta. Inverse of the temperature vector in parallel tempering
%       - opt. The PRML options
%
% Return:
%       - states. Same
%------------------------------------------------------------------------------


            a = ceil(rand*opt.Nchain);

            if a == opt.Nchain
                b = 1;
            else
                b = a+1;
            end

            SSa = states(a, opt.Nparams+1);
            SSb = states(b, opt.Nparams+1);

            % Chains are swaped with certain Metropolis probability
            rho = (exp(-0.5*(SSa-SSb)/sigmaSqu(a)))^(Beta(b)-Beta(a));

            if rand < min(1, rho)
                temp         = states(a, :);
                states(a, :) = states(b, :);
                states(b, :) = temp;
                clear temp;
            end

        end

        function Temperatures = temperatureDynamics(it, states, Temperatures, sigmaSqu, opt)
%------------------------------------------------------------------------------
% Temperature evolution in the parallel tempering
%
% parameter:
%       - it. The index of the current chain
%       - states. The Nchain x (Nparams+1) matrix
%       - Temperatures. The vector of the temperature in the parallel tempering
%       - simgaSqu. The sigma square vector (the covariance matrix)
%       - opt. The PRML options
%
% Return:
%       - Temperatures. Same
%------------------------------------------------------------------------------


            t0 = 1e3; nu = 100;

            Beta = 1 ./ Temperatures(it, :);

            for i = 2: opt.Nchain

                b = i - 1;
                if i == opt.Nchain
                    c = 1;
                else
                    c = i + 1;
                end

                SSa = states(i, opt.Nparams+1);
                SSb = states(b, opt.Nparams+1);
                SSc = states(c, opt.Nparams+1);

                rho_ab = min(1, (exp(-0.5*(SSa-SSb)/sigmaSqu(i)))^(Beta(b)-Beta(i)) );
                rho_ca = min(1, (exp(-0.5*(SSc-SSa)/sigmaSqu(c)))^(Beta(i)-Beta(c)) );

                differential = t0/(nu*(it+t0)) * (rho_ab - rho_ca);

                Temperatures(it+1, i) = Temperatures(it+1, i-1) + ...
                    exp(log(Temperatures(it, i)-Temperatures(it, i-1)) + differential);

            end

        end

        function Population = convertASCIItoMAT_PRML(maxIter, opt)
%------------------------------------------------------------------------------
% Load and convert the data in ascii format to the mat format
%   then generate the population for statistics
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgorithms.convertASCIItoMAT_PRML: There are no enough input arguments \n');
            end

            chain = [];

            for i = 1:opt.Nchain

                load(sprintf('chainData_%d.dat', i));

                chain = [chain eval(sprintf('chainData_%d', i))];

            end

            save('chain.dat', 'chain', '-ascii');
            for  j = 1:opt.Nchain
                eval(sprintf('delete chainData_%d.dat',j));
            end

            % Discard the former 50% chain, retain only the last 50%
            if maxIter < opt.nsamples
                idx = floor(0.5 * maxIter);
            else
                idx = floor(0.5 * opt.nsamples);
            end

            for k = 1:1
                eval(sprintf('chainData_%d(1:idx, :) = [];', k));
                Population = [Population; eval(sprintf('chainData_%d', k))];
            end

            save('population.dat', 'Population', '-ascii');

        end


    end %PRML


%   MADE
    methods (Static = true, Access = 'public')

        function [xValue, yValue] = Metropolis_Adjusted_Differential_Evolution(obj, params)
% -----------------------------------------------------------------------------
% Metropolis Adjusted Differential Evolution algorithm (MADE)
%
% MADE optimizes a problem by combining the prominent features of Metropolis
% Hastings algorithm and Differential Evolution algorithm. In the upper level,
% each chain is accepted with the Metropolis probability, while in the lower
% level, chains have an evolution with resort to heuristic method, Differential
% Evolution algorithm.
%
% Unlike the algorithm, PSO and DE, the MADE obtain the parameter distributions
% rather than the single parameter set. It provides the confidential intervals
% for each parameter, since it is based on the Markov Chain Monte Carlo (MCMC).
%
%  Parameter:
%        - params. It is the specified parameters from the main function.
%        And let the main function informed that which parameters are need
%        to be optimized
% -----------------------------------------------------------------------------


            startTime = clock;

            if nargin < 2
                error('OptAlgorithms.MADE: There are no enough input arguments \n');
            end

            % Get the sampler options
            opt = OptAlgorithms.getOptions_MADE(obj, params);

            % Initialization of the chains
            states = OptAlgorithms.initChains(opt);

            % Preallocation
            accepted = 0; n0 = 1;

            chain = zeros(opt.Nchain, opt.Nparams+1, opt.nsamples);

            % Calculate the initial sigma square values
            sumSquare = states(:, opt.Nparams+1);
            if any(sumSquare < 0)
                row = find(sumSquare < 0);
                sumSquare(row) = exp( sumSquare(row) );
            end
            sigmaSqu = sumSquare ./ (opt.nDataPoint - opt.Nparams);

            sigmaSqu_0 = sigmaSqu;

%-----------------------------------------------------------------------------------------

% The main loop
%-----------------------------------------------------------------------------------------
            for i = 1: opt.nsamples

                fprintf('Iter: %4d ----- accept_ratio: %3d%%', i, fix(accepted/(opt.Nchain*i)*100));

                [states, accepted] = OptAlgorithms.MADE_sampler(states, sigmaSqu, accepted, opt);

                % Append the newest states to the end of the 3-D matrix
                chain(:,:,i) = states;

                for j = 1:opt.Nchain
                    temp = states(j, :);
                    save(sprintf('chainData_%d.dat', j), 'temp', '-ascii', '-append');
                end

                % In each opt.convergInt interval, check the convergence condition
                if mod(i, opt.convergInt) == 0 && i > opt.convergInt

                    criterion = OptAlgorithms.GelmanR_statistic(i, chain(:,:,1:i), opt);

                    if all(criterion < opt.criterionTol) || i == opt.nsamples
                        maxIter = i;
                        break
                    end

                end

                % Variance of error distribution (sigma) was treated as a parameter to be estimated.
                sumSquare = states(:, opt.Nparams+1);
                if any(sumSquare < 0)
                    row = find(sumSquare < 0);
                    sumSquare(row) = exp( sumSquare(row) );
                end

                for k = 1:opt.Nchain
                    sigmaSqu(k)  = 1 ./ OptAlgorithms.GammarDistribution(1, 1, (n0 + opt.nDataPoint)/2, ...
                        2 / (n0 * sigmaSqu_0(k) + sumSquare(k)));
                end

            end
%-----------------------------------------------------------------------------------------

% Post-process
%-----------------------------------------------------------------------------------------

            clear chain;

            % Generate the population for figure plot
            maxIter = 1000
            Population = OptAlgorithms.convertASCIItoMAT_MADE(maxIter, opt);

            % Plot the population
            OptAlgorithms.FigurePlot(Population, opt);

            [yValue, row]  = min(Population(:,opt.Nparams+1));
            xValue         = Population(row,1:opt.Nparams);

            % Gather some useful information and store them
            result.optTime        = etime(clock,startTime) / 3600;
            result.convergDiagno  = criterion;
            result.Nchain         = opt.Nchain;
            result.Ieteration     = maxIter;
            result.acceptanceRate = fix( accepted / (opt.Nchain * maxIter) * 100 );
            result.correlation    = corrcoef(Population(:,1:opt.Nparams));
            result.population     = Population;
            result.xValue         = xValue;
            result.yValue         = yValue;

            save(sprintf('result_%3d.mat', fix(rand*1000)), 'result');
            fprintf('Markov chain simulation finished and the result saved as result.mat \n');

        end

        function opt = getOptions_MADE(obj, params)
% -----------------------------------------------------------------------------
%  The parameters for the Optimizer
%
%  Parameter:
%       - params. It is the specified parameters from the main function.
%        And let the main function informed that which parameters are need
%        to be optimized.
%
%  Return:
%       - opt.
%           + Nchain. The number of the candidates (particles)
%           + Nparams. The number of the optimized parameters
%           + bounds. The boundary limitation of the parameters
%           + nsamples. The maximum number of algorithm's iteration
%           + cr_Probability. The cross-over probability in DE's formula
%           + weight. The weight coefficients in DE's formula
%           + strategy. There are 1,2,3,4,5,6 strategies in DE algorithm to 
%               deal with the cross-over and mutation. As for detais, we
%               will refer your the original paper of Storn and Price
% -----------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgorithms.getOptions_MADE: There are no enough input arguments \n');
            end

            opt = [];

            opt.Nchain        = OptAlgorithms.swarm;
            opt.Nparams       = length(fieldnames(params));
            opt.bounds        = obj.paramBound';

            opt.nsamples      = OptAlgorithms.sample;
            opt.criterionTol  = 1.01;
            opt.convergInt    = 50;
            % this data should be manually changed, as it is related to the data
            opt.nDataPoint    = OptAlgorithms.dataPoint;

            % Check out the dimension of the set of parameters, and the boundary limitation
            [row, col] = size(opt.Nparams);
            if row > 1 || col > 1
                error('OptAlgorithms.getOptions_MADE: The initialized dimension of the set of parameters might be wrong \n');
            end

            [row, col] = size(opt.bounds);
            if col ~= opt.Nparams || row ~= 2
                error('OptAlgorithms.getOptions_MADE: Please check your setup of the range of parameters \n');
            end

            % Those are the specific parameters for the DE kernel
            opt.cr_Probability  = 0.5;
            opt.weight          = 0.3;
            opt.strategy        = OptAlgorithms.DE_strategy;

        end

        function initChain = initChains(opt)
% -----------------------------------------------------------------------------
% The initilization of the chains
%
%  Parameter:
%       - opt.
%           + PopulSize. The number of the candidates (particles)
%           + IndivSize. The number of the optimized parameters
%           + IndivScope. The boundary limitation of the parameters
%           + loopCount. The maximum number of algorithm's iteration
%           + cr_Probability. The cross-over probability in DE's formula
%           + weight. The weight coefficients in DE's formula
%           + strategy. There are 1,2,3,4,5,6 strategies in DE algorithm to 
%               deal with the cross-over and mutation. As for detais, we
%               will refer your the original paper of Storn and Price
%
% Return:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
% -----------------------------------------------------------------------------


        %   Initilization of the chains
            initChain = rand(opt.Nchain, opt.Nparams+1);

        %   Use vectorization to speed up. Keep the parameters in the domain
            initChain(:,1:opt.Nparams) = repmat(opt.bounds(1,:),opt.Nchain,1) + ...
                initChain(:,1:opt.Nparams) .* repmat( (opt.bounds(2,:) - opt.bounds(1,:)),opt.Nchain,1 );

        %   Simulation of the sampled points
            initChain(:,opt.Nparams+1) = arrayfun( @(idx) feval( OptAlgorithms.FUNC, ...
                initChain(idx,1:opt.Nparams) ), 1: opt.Nchain);

        end

        function [states, accepted] = MADE_sampler(states, sigmaSqu, accepted, opt)
%-----------------------------------------------------------------------------------------
% The DE sampler for the Metropolis adjusted differential evolution method
%-----------------------------------------------------------------------------------------


            if nargin < 4
                error('OptAlgorithms.MADE_sampler: There are no enough input arguments \n');
            end

            % The evolution of the chains with DE kernel
            [tempPopulation, OptPopul] = OptAlgorithms.MADE_DE_Evolution(states, opt);

            % Abstract best information so far from the population and display it
            fprintf('----------------  Minimum: %g  ---------------- \n', OptPopul(opt.Nchain+1, opt.Nparams+1));
            fprintf('%10.3g | ', OptPopul(1, 1:opt.Nparams)); fprintf('\n');

            % In each chain, the proposal point is accepted in terms of the Metropolis probability
            for j = 1: opt.Nchain

                SS = states(j, opt.Nparams+1);

                proposal = tempPopulation(j, 1:opt.Nparams);

                if any(proposal < opt.bounds(1,:)) || any(proposal > opt.bounds(2,:))
                    newSS = Inf;
                else
                    newSS = feval( OptAlgorithms.FUNC, proposal );
                end

                rho = exp( -0.5*(newSS - SS) / sigmaSqu(j));

                if rand <= min(1, rho)
                    states(j, 1:opt.Nparams) = proposal;
                    states(j, opt.Nparams+1) = newSS;
                    accepted = accepted + 1;
                end

            end

        end

        function [tempPop, OptPopul] = MADE_DE_Evolution(Population, opt)
% -----------------------------------------------------------------------------
% The evolution of population
%
% Parameters:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
%       - opt. Please see the comments of the function, InitPopulation
%
% Return:
%       - Population. The population of the particles, correspondingly the objective value
%       - OptPopul. The best fit found among the population.
% -----------------------------------------------------------------------------


            R = opt.Nchain;
            C = opt.Nparams;
            [minValue, row] = min(Population(:,opt.Nparams+1));

            OptPopul = zeros(R+1, C+1);
            OptPopul(R+1, C+1) = minValue;

            OptPopul(1:R, 1:C) = repmat( Population(row, 1:C), R, 1 );
            clear row;

            indexArray = (0:1:R-1);
            index = randperm(4);

            cr_mutation = rand(R, C) < opt.cr_Probability;
            cr_old = cr_mutation < 0.5;

            ShuffRow1 = randperm(R);
            PopMutR1 = Population(ShuffRow1, 1:C);

            idxShuff = rem(indexArray + index(1), R);
            ShuffRow2 = ShuffRow1(idxShuff + 1);
            PopMutR2 = Population(ShuffRow2, 1:C);

            idxShuff = rem(indexArray + index(2), R);
            ShuffRow3 = ShuffRow2(idxShuff + 1);
            PopMutR3 = Population(ShuffRow3, 1:C);

            idxShuff = rem(indexArray + index(3), R);
            ShuffRow4 = ShuffRow3(idxShuff + 1);
            PopMutR4 = Population(ShuffRow4, 1:C);

            idxShuff = rem(indexArray + index(4), R);
            ShuffRow5 = ShuffRow4(idxShuff + 1);
            PopMutR5 = Population(ShuffRow5, 1:C);

            switch opt.strategy
                case 1
                % strategy 1
                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * opt.weight;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 2
                % strategy 2
                    tempPop = Population(1:R, 1:C) + opt.weight * (OptPopul(1:R, 1:C) - Population(1:R, 1:C)) + (PopMutR1 - PopMutR2) * opt.weight;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 3
                % strategy 3
                    tempPop = OptPopul(1:R, 1:C) + (PopMutR1 - PopMutR2) .* ((1 -0.9999) * rand(R, C) + opt.weight);
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 4
                % strategy 4
                    f1 = (1-opt.weight) * rand(R, 1) + opt.weight;
                    for i = 1: C
                        PopMutR5(:, i) = f1;
                    end

                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) .* PopMutR5;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 5
                % strategy 5
                    f1 = (1-opt.weight) * rand + opt.weight;
                    tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * f1;
                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;

                case 6
                % strategy 6
                    if (rand < 0.5)
                        tempPop = PopMutR3 + (PopMutR1 - PopMutR2) * opt.weight;
                    else
                        tempPop = PopMutR3 + 0.5 * (opt.weight + 1.0) * (PopMutR1 + PopMutR2 - 2 * PopMutR3);
                    end

                    tempPop = Population(1:R, 1:C) .* cr_old + tempPop .* cr_mutation;
                 % add new
            end


            % Check the boundary limitation
            loBound = repmat(opt.bounds(1,:), R, 1);
            [row, col] = find( (tempPop(1:R, 1:C) - loBound) < 0 );
            tempPop((col-1).*R + row) = loBound((col-1).*R + row);

            upBound = repmat(opt.bounds(2,:), R, 1);
            [row, col] = find( (tempPop(1:R, 1:C) - upBound) > 0 );
            tempPop((col-1).*R + row) = upBound((col-1).*R + row);
            clear row col;


        end  

        function y = GammarDistribution(m, n, a, b)
%-----------------------------------------------------------------------------------------
% GammarDistrib random deviates from gamma distribution
% 
%  GAMMAR_MT(M,N,A,B) returns a M*N matrix of random deviates from the Gamma
%  distribution with shape parameter A and scale parameter B:
%
%  p(x|A,B) = B^-A/gamma(A)*x^(A-1)*exp(-x/B)
%
%  Uses method of Marsaglia and Tsang (2000)
%
% G. Marsaglia and W. W. Tsang:
% A Simple Method for Generating Gamma Variables,
% ACM Transactions on Mathematical Software, Vol. 26, No. 3, September 2000, 363-372.
%-----------------------------------------------------------------------------------------

            if nargin < 4, b = 1; end

            y = zeros(m, n);
            for j=1: n
                for i=1: m
                    y(i, j) = OptAlgorithms.Gammar(a, b);
                end
            end

        end

        function y = Gammar(a, b)

            if a < 1

                y = OptAlgorithms.Gammar(1+a, b) * rand(1) ^ (1/a);

            else

                d = a - 1/3;
                c = 1 / sqrt(9*d);

                while(1)

                    while(1)
                        x = randn(1);
                        v = 1 + c*x;

                        if v > 0, break, end

                    end

                    v = v^3; u = rand(1);

                    if u < 1 - 0.0331*x^4, break, end

                    if log(u) < 0.5 * x^2 + d * (1-v+log(v)), break, end

                end

                y = b * d * v;

            end

        end

        function criterion = GelmanR_statistic(idx, chain, opt)
%-----------------------------------------------------------------------------------------
%
%-----------------------------------------------------------------------------------------


            if nargin < 3
                error('OptAlgorithms.GelmanR_statistic: There are no enough input arguments \n');
            end

            % Split each chain into half and check all the resulting half-sequences
            index           = floor(0.5 * idx); 
            eachChain       = zeros(index, opt.Nparams);
            betweenMean     = zeros(opt.Nchain, opt.Nparams);
            withinVariance  = zeros(opt.Nchain, opt.Nparams);

            % Mean and variance of each half-sequence chain
            for i = 1: opt.Nchain

                for j = 1: opt.Nparams
                    for k = 1: index
                        eachChain(k,j) = chain(i,j,k+index);
                    end
                end

                betweenMean(i,:)    = mean(eachChain);
                withinVariance(i,:) = var(eachChain);

            end

            % Between-sequence variance
            Sum = 0;
            for i = 1: opt.Nchain
               Sum = Sum + (betweenMean(i,:) - mean(betweenMean)) .^ 2;
            end
            B = Sum ./ (opt.Nchain-1);

            % Within-sequence variance
            Sum = 0;
            for i = 1: opt.Nchain
                Sum = Sum + withinVariance(i,:);
            end
            W = Sum ./ opt.Nchain;

            % Convergence diagnostics
            criterion = sqrt(1 + B ./ W);

        end

        function Population = convertASCIItoMAT_MADE(maxIter, opt)
%------------------------------------------------------------------------------
% Load and convert the data in ascii format to the mat format
%   then generate the population for statistics
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgorithms.convertASCIItoMAT_MADE: There are no enough input arguments \n');
            end

            chain = []; Population = [];

            for i = 1:opt.Nchain

                load(sprintf('chainData_%d.dat', i));

                chain = [chain eval(sprintf('chainData_%d', i))];

            end

            save('chain.dat', 'chain', '-ascii');
            for  j = 1:opt.Nchain
                eval(sprintf('delete chainData_%d.dat',j));
            end

            % Discard the former 50% chain, retain only the last 50%
            if maxIter < opt.nsamples
                idx = floor(0.5 * maxIter);
            else
                idx = floor(0.5 * opt.nsamples);
            end

            for k = 1:opt.Nchain
                eval(sprintf('chainData_%d(1:idx, :) = [];', k));
                Population = [Population; eval(sprintf('chainData_%d', k))];
            end

            save('population.dat', 'Population', '-ascii', '-append');

        end

        function FigurePlot(Population, opt)
%------------------------------------------------------------------------------
% Plot the histgram and scatter figures of the last population
%------------------------------------------------------------------------------


            if nargin < 2
                error('OptAlgorithms.FigurePlot: There are no enough input arguments \n');
            end

            figure(1);clf
            for i = 1: opt.Nparams

                subplot(round(opt.Nparams/2),2,i,'Parent', figure(1));

                histfit( Population(:,i), 50, 'kernel');

                xlabel(sprintf('$x_%d$', i), 'FontSize', 20, 'Interpreter', 'latex');
                ylabel(sprintf('Frequency'), 'FontSize', 20, 'Interpreter', 'latex');
                set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
%                OptAlgorithms.tickLabelFormat(gca, 'x', '%0.2e');
%                OptAlgorithms.tickLabelFormat(gca, 'x', []);
                set(gca, 'XTickLabel', num2str(get(gca, 'xTick')', '%3g'));
%                 OptAlgorithms.xtickLabelRotate([], 15, [], 'FontSize', 20, 'FontName', 'Times New Roman');
                set(gca, 'ygrid', 'on');

            end

            figure(2);clf
            for i = 1: opt.Nparams-1
                for j = i: opt.Nparams-1

                    subplot(opt.Nparams-1, opt.Nparams-1, j+(i-1)*(opt.Nparams-1), 'Parent', figure(2));

                    scatter(Population(:,j+1), Population(:,i));

                    xlabel(sprintf('$x_%d$', j+1), 'FontName', 'Times New Roman', 'FontSize', 20, 'Interpreter', 'latex');
                    ylabel(sprintf('$x_%d$', i), 'FontName', 'Times New Roman', 'FontSize', 20, 'Interpreter', 'latex');
                    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
%                    OptAlgorithms.tickLabelFormat(gca, 'x', '%0.2e');
                    set(gca, 'XTickLabel', num2str(get(gca, 'xTick')', '%3g'));
%                     OptAlgorithms.xtickLabelRotate([], 15, [], 'FontSize', 20, 'FontName', 'Times New Roman');
                    grid on;

                end
           end

        end

    end % MADE


    methods (Static = true, Access = 'public')

        function checkOptDimension(opt, LEN)
%------------------------------------------------------------------------------
% Check the consistency of the parameter set and the boundary set
%------------------------------------------------------------------------------


            if length(opt.paramBound) ~= LEN
                error('The setup of the initial boundary in SMBOptimiztion is incorrent \n');
            end
        end

        function tickLabelFormat(hAxes, axName, format)
%------------------------------------------------------------------------------
% Sets the format of the tick labels
%
% Syntax:
%    ticklabelformat(hAxes,axName,format)
%
% Input Parameters:
%    hAxes  - handle to the modified axes, such as returned by the gca function
%    axName - name(s) of axles to modify: 'x','y','z' or combination (e.g. 'xy')
%    format - format of the tick labels in sprintf format (e.g. '%.1f V') or a
%             function handle that will be called whenever labels need to be updated
%
%    Note: Calling TICKLABELFORMAT again with an empty ([] or '') format will revert
%    ^^^^  to Matlab's normal tick labels display behavior
%
% Examples:
%    ticklabelformat(gca,'y','%.6g V') - sets y axis on current axes to display 6 significant digits
%    ticklabelformat(gca,'xy','%.2f')  - sets x & y axes on current axes to display 2 decimal digits
%    ticklabelformat(gca,'z',@myCbFcn) - sets a function to update the Z tick labels on current axes
%    ticklabelformat(gca,'z',{@myCbFcn,extraData}) - sets an update function as above, with extra data
%
% Warning:
%    This code heavily relies on undocumented and unsupported Matlab functionality.
%    It works on Matlab 7+, but use at your own risk!
%
% Technical description and more details:
%    http://UndocumentedMatlab.com/blog/setting-axes-tick-labels-format
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany@gmail.com)
%------------------------------------------------------------------------------


            % Check # of args (we now have narginchk but this is not available on older Matlab releases)
            if nargin < 1
                help(mfilename)
                return;
            elseif nargin < 3
                error('OptAlgorithms.ticklabelformat: Not enough input arguments \n');
            end

            % Check input args
            if ~ishandle(hAxes) || (~isa(handle(hAxes),'axes') && ~isa(handle(hAxes),'matlab.graphics.axis.Axes'))
                error('OptAlgorithms.ticklabelformat: hAxes input argument must be a valid axes handle \n');
            elseif ~ischar(axName)
                error('OptAlgorithms.ticklabelformat: axName input argument must be a string \n');
            elseif ~isempty(format) && ~ischar(format) && ~isa(format,'function_handle') && ~iscell(format)
                error('OptAlgorithms.ticklabelformat: format input argument must be a string or function handle \n');
            end

            % normalize axes name(s) to lowercase
            axName = lower(axName);

            if strfind(axName,'x')
                install_adjust_ticklbl(hAxes,'X',format)
            elseif strfind(axName,'y')
                install_adjust_ticklbl(hAxes,'Y',format)
            elseif strfind(axName,'z')
                install_adjust_ticklbl(hAxes,'Z',format)
            end

            function install_adjust_ticklbl(hAxes,axName,format)
            % Install the new tick labels for the specified axes

                % If empty format was specified
                if isempty(format)

                    % Remove the current format (revert to default Matlab format)
                    set(hAxes,[axName 'TickLabelMode'],'auto')
                    setappdata(hAxes,[axName 'TickListener'],[])
                    return

                end

                % Determine whether to use the specified format as a
                % sprintf format or a user-specified callback function
                if ischar(format)
                    cb = {@adjust_ticklbl axName format};
                else
                    cb = format;
                end

                % Now install axis tick listeners to adjust tick labels
                % (use undocumented feature for adjustments)
                ha = handle(hAxes);
                propName = [axName 'Tick'];
                hp = findprop(ha,propName);

                try
                    % R2014a or older
                    hl = handle.listener(ha,hp,'PropertyPostSet',cb);

                    % Adjust tick labels now
                    % eventData.AffectedObject = hAxes;
                    % adjust_ticklbl([],eventData,axName,format)
                    set(hAxes,[axName 'TickLabelMode'],'manual')
                    set(hAxes,[axName 'TickLabelMode'],'auto')
                catch
                    % R2014b or newer
                    if iscell(cb)
                        cb = @(h,e) feval(cb{1},h,e,cb{2:end});
                    end
                    % hl(1) = addlistener(ha,propName,'PostSet',cb);
                    hl(1) = event.proplistener(ha,hp,'PostSet',cb);

                    % *Tick properties don't trigger PostSet events when updated automatically in R2014b - need to use *Lim
                    % addlistener(ha,[axName 'Lim'],'PostSet',cb);
                    hRuler = get(ha,[axName 'Ruler']);
                    % hl(2) = addlistener(hRuler,'MarkedClean',cb);
                    hl(2) = event.listener(hRuler,'MarkedClean',cb);

                    % Adjust tick labels now
                    eventData.AffectedObject = hAxes;
                    if ischar(format)
                        adjust_ticklbl([],eventData,axName,format)
                    else
                        hgfeval(format);
                    end
                    set(hAxes,[axName 'TickLabelMode'],'manual')
                    % set(hAxes,[axName 'TickLabelMode'],'auto')  % causes labels not to be updated in R2014b!
                end

                setappdata(hAxes,[axName 'TickListener'],hl)
                % drawnow;

                function adjust_ticklbl(hProp,eventData,axName,format)
                % Default tick labels update callback function (used if user did not specify their own function)

                    try
                        hAxes = eventData.AffectedObject;
                    catch
                        hAxes = ancestor(eventData.Source,'Axes');
                    end
                    tickValues = get(hAxes,[axName 'Tick']);
                    tickLabels = arrayfun(@(x)(sprintf(format,x)),tickValues,'UniformOutput',false);
                    set(hAxes,[axName 'TickLabel'],tickLabels);

                end

            end

        end

        function hText = xtickLabelRotate(XTick, rot, varargin)
%------------------------------------------------------------------------------
% xtick label rotate
%
% Parameter:
%       - XTick. vector array of XTick positions & values (numeric) uses current 
%           XTick values or XTickLabel cell array by default (if empty) 
%       - rot. angle of rotation in degrees, 90° by default
%       - XTickLabel: cell array of label strings
%       - [var]. Optional. "Property-value" pairs passed to text generator
%           ex: 'interpreter','none', 'Color','m','Fontweight','bold'
%
% Return:
%       - hText. handle vector to text labels
%
% Example 1:  Rotate existing XTickLabels at their current position by 90°
%    xticklabel_rotate
%
% Example 2:  Rotate existing XTickLabels at their current position by 45° and change
% font size
%    xticklabel_rotate([],45,[],'FontSize',14)
%
% Example 3:  Set the positions of the XTicks and rotate them 90°
%    figure;  plot([1960:2004],randn(45,1)); xlim([1960 2004]);
%    xticklabel_rotate([1960:2:2004]);
%
% Example 4:  Use text labels at XTick positions rotated 45° without tex interpreter
%    xticklabel_rotate(XTick,45,NameFields,'interpreter','none');
%
% Example 5:  Use text labels rotated 90° at current positions
%    xticklabel_rotate([],90,NameFields);
%
% Example 6:  Multiline labels
%    figure;plot([1:4],[1:4])
%    axis([0.5 4.5 1 4])
%    xticklabel_rotate([1:4],45,{{'aaa' 'AA'};{'bbb' 'AA'};{'ccc' 'BB'};{'ddd' 'BB'}})
%
% This is a modified version of xticklabel_rotate90 by Denis Gilbert
% Modifications include Text labels (in the form of cell array)
%       Arbitrary angle rotation
%       Output of text handles
%       Resizing of axes and title/xlabel/ylabel positions to maintain same overall size 
%          and keep text on plot
%          (handles small window resizing after, but not well due to proportional placement with 
%           fixed font size. To fix this would require a serious resize function)
%       Uses current XTick by default
%       Uses current XTickLabel is different from XTick values (meaning has been already defined)
%
% Author: Brian FG Katz, bfgkatz@hotmail.com
%------------------------------------------------------------------------------

            % check to see if xticklabel_rotate has already been here (no other reason for this to happen)
            if isempty(get(gca,'XTickLabel')),
                error('OptAlgorithms.xtickLabelRotate: can not process, either xticklabel_rotate has already been run or XTickLabel field has been erased \n');
            end

            % Modified with forum comment by "Nathan Pust" allow the current text labels to be used and property value pairs to be changed for those labels
            if (nargin < 3 || isempty(varargin{1})) && (~exist('XTick') || isempty(XTick))

                xTickLabels = get(gca,'XTickLabel');

                if ~iscell(xTickLabels)
                    % remove trailing spaces if exist (typical with auto generated XTickLabel)
                    temp1 = num2cell(xTickLabels,2);
                    for loop = 1:length(temp1),
                        temp1{loop} = deblank(temp1{loop});
                    end
                    xTickLabels = temp1;
                end

                varargin = varargin(2:length(varargin));

            end

            % if no XTick is defined use the current XTick
            if (~exist('XTick') || isempty(XTick))
                XTick = get(gca,'XTick');
            end

		    % Make XTick a column vector
            XTick = XTick(:);

            if ~exist('xTickLabels')
                % Define the xtickLabels 
                % If XtickLabel is passed as a cell array then use the text
                if ~isempty(varargin) && (iscell(varargin{1}))
                    xTickLabels = varargin{1};
                    varargin = varargin(2:length(varargin));
                else
                    xTickLabels = num2str(XTick);
                end
            end

            if length(XTick) ~= length(xTickLabels)
                error('OptAlgorithms.xtickLabelRotate: must have same number of elements in "XTick" and "XTickLabel" \n');
            end

            % Set the Xtick locations and set XTicklabel to an empty string
            set(gca,'XTick',XTick,'XTickLabel','');

            if nargin < 2
                rot = 90 ;
            end

            % Determine the location of the labels based on the position of the xlabel
            hxLabel = get(gca,'XLabel');
            xLabelString = get(hxLabel,'String');

            set(hxLabel,'Units','data');
            xLabelPosition = get(hxLabel,'Position');
            y = xLabelPosition(2);

            % CODE below was modified following suggestions from Urs Schwarz
            y = repmat(y,size(XTick,1),1);
            % retrieve current axis' fontsize
            fs = get(gca,'fontsize');

            if ~iscell(xTickLabels)
                % Place the new xTickLabels by creating TEXT objects
                hText = text(XTick, y, xTickLabels,'fontsize',fs);
            else
                % Place multi-line text approximately where tick labels belong
                for cnt=1:length(XTick)
                    hText(cnt) = text(XTick(cnt),y(cnt),xTickLabels{cnt}, ...
                        'VerticalAlignment','top', 'UserData','xtick');
                end
            end

            % Rotate the text objects by ROT degrees
            % Modified with modified forum comment by "Korey Y" to deal with labels at top
            % Further edits added for axis position
            xAxisLocation = get(gca, 'XAxisLocation');
            if strcmp(xAxisLocation,'bottom')
                set(hText,'Rotation',rot,'HorizontalAlignment','right',varargin{:});
            else
                set(hText,'Rotation',rot,'HorizontalAlignment','left',varargin{:});
            end

            % Adjust the size of the axis to accomodate for longest label (like if they are text ones)
            % This approach keeps the top of the graph at the same place and tries to keep xlabel at the same place
            % This approach keeps the right side of the graph at the same place 

            set(get(gca,'xlabel'),'units','data');
                labxorigpos_data = get(get(gca,'xlabel'),'position');
            set(get(gca,'ylabel'),'units','data');
                labyorigpos_data = get(get(gca,'ylabel'),'position');
            set(get(gca,'title'),'units','data');
                labtorigpos_data = get(get(gca,'title'),'position');

            set(gca,'units','pixel');
            set(hText,'units','pixel');
            set(get(gca,'xlabel'),'units','pixel');
            set(get(gca,'ylabel'),'units','pixel');

            origpos = get(gca,'position');

            % Modified with forum comment from "Peter Pan" to deal with case when only one XTickLabelName is given. 
            x = get( hText, 'extent' );
            if iscell( x ) == true
                textsizes = cell2mat( x );
            else
                textsizes = x;
            end

            largest =  max(textsizes(:,3));
            longest =  max(textsizes(:,4));

            laborigext = get(get(gca,'xlabel'),'extent');
            laborigpos = get(get(gca,'xlabel'),'position');

            labyorigext = get(get(gca,'ylabel'),'extent');
            labyorigpos = get(get(gca,'ylabel'),'position');
            leftlabdist = labyorigpos(1) + labyorigext(1);

            % assume first entry is the farthest left
            leftpos = get(hText(1),'position');
            leftext = get(hText(1),'extent');
            leftdist = leftpos(1) + leftext(1);
            if leftdist > 0, leftdist = 0; end

            % Modified to allow for top axis labels and to minimize axis resizing
            if strcmp(xAxisLocation,'bottom')
                newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
                    origpos(2)+((longest+laborigpos(2))-get(gca,'FontSize')) ...
                    origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
                    origpos(4)-((longest+laborigpos(2))-get(gca,'FontSize'))];
            else
                newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
                    origpos(2) ...
                    origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
                    origpos(4)-(longest)+get(gca,'FontSize')];
            end
            set(gca,'position',newpos);

            % readjust position of text labels after resize of plot
            set(hText,'units','data');
            for loop= 1:length(hText)
                set(hText(loop),'position',[XTick(loop), y(loop)]);
            end

            % adjust position of xlabel and ylabel
            laborigpos = get(get(gca,'xlabel'),'position');
            set(get(gca,'xlabel'),'position',[laborigpos(1) laborigpos(2)-longest 0]);

            % switch to data coord and fix it all
            set(get(gca,'ylabel'),'units','data');
            set(get(gca,'ylabel'),'position',labyorigpos_data);
            set(get(gca,'title'),'position',labtorigpos_data);

            set(get(gca,'xlabel'),'units','data');
                labxorigpos_data_new = get(get(gca,'xlabel'),'position');
            set(get(gca,'xlabel'),'position',[labxorigpos_data(1) labxorigpos_data_new(2)]);

            % Reset all units to normalized to allow future resizing
            set(get(gca,'xlabel'),'units','normalized');
            set(get(gca,'ylabel'),'units','normalized');
            set(get(gca,'title'),'units','normalized');
            set(hText,'units','normalized');
            set(gca,'units','normalized');

            if nargout < 1,
                clear hText
            end

        end

    end


% Upper level algorithm, in charge of discrete structural optimization
% -----------------------------------------------------------------------------
    methods (Static = true, Access = 'public')

        function structure = discreteInit(opt)
% -----------------------------------------------------------------------------
% Generating the initial population of structures
%
% Encoding: each node between two adjacent columns are denoted by a number sequently
%           0 1 2 3 4 5 6 7 8 9 10 ...
% Here 0 represents the desorbent port all the time. As it is a loop, we always need a starting point
% And the sequence of ports are D E (ext_1 ext_2) F R constantly. The selective ranges of each pointer
% (E,F,R) are shown as follows in the binary scenario:
%           0 1 2 3 4 5 6 7 8 9 10 ...
%           D 
%             E < ------- > E      : extract_pool
%               F < ------- > F    : feed_pool
%                 R < ------- > R  : raffinate_pool
% -----------------------------------------------------------------------------


            nodeIndex = opt.nColumn -1 ;

            % Preallocate of the structure matrix
            structure = zeros(opt.structNumber,opt.nZone+1);

            for i = 1:opt.structNumber

                if opt.nZone == 4
                    extract_pool = [1, nodeIndex-2];
                    structure(i,2) = randi(extract_pool);

                    feed_pool = [structure(i,2)+1, nodeIndex-1];
                    structure(i,3) = randi(feed_pool);

                    raffinate_pool = [structure(i,3)+1, nodeIndex];
                    structure(i,4) = randi(raffinate_pool);

                elseif opt.nZone == 5
                    extract1_pool = [1, nodeIndex-3];
                    structure(i,2) = randi(extract1_pool);

                    extract2_pool = [structure(i,2)+1, nodeIndex-2];
                    structure(i,3) = randi(extract2_pool);

                    feed_pool = [structure(i,3)+1, nodeIndex-1];
                    structure(i,4) = randi(feed_pool);

                    raffinate_pool = [structure(i,4)+1, nodeIndex];
                    structure(i,5) = randi(raffinate_pool);

                end
            end

        end

        function structID = structure2structID(opt,structure)
% -----------------------------------------------------------------------------
% This is the rountine that decode the structure into structID for simulation
% 
% For instance, the structure is [0, 3, 5, 9] in a binary situation with column amount 10
% the structID is [3,2,4,1], there are three in the zone I, two in zone II, four in
% zone III, one in zone IV
% -----------------------------------------------------------------------------


            structID = zeros(1, opt.nZone);

            structID(1:opt.nZone-1) = structure(2:end) - structure(1:end-1);

            structID(end) = opt.nColumn - structure(end);

        end

        function [theta, objective] = continuousUnitOptimization(opt, params, optimization_method)
% -----------------------------------------------------------------------------
% This is the main function of the optimization of the Simulated Moving Bed
% The optimized parameters in this case are
%       - columnLength
%       - switchTime
%       - flowRates_recycle
%       - flowRate_feed
%       - flowRate_desorbent
%       - flowRate_extract
%
%       theta = {L_c, t_s, Q_{re}, Q_F, Q_D, Q_E}
%
% In the FIVE-ZONE, the optimized parameters are
%       - columnLength
%       - switchTime
%       - flowRates_recycle
%       - flowRate_feed
%       - flowRate_desorbent
%       - flowRate_extract_1
%       - flowRate_extract_2
%
%       theta = {L_c, t_s, Q_{re}, Q_F, Q_D, Q_{E1}, Q_{E2}}
%
% There are four types of algorithms are integrated into this code, either
% based on Heuristical theory or Deterministic theory, either optimization or sampling.
%       - Particle Swarm Optimizatio (PSO)
%       - Differential Evolution (DE)
%       - Metropolis Adjusted Differential Evolution (MADE)
%       - Markov Chain Monte Carlo (MCMC)
%       - Metropolis Adjusted Langevin Algorithm (MLA) defined on Riemann manifold
% 
% -----------------------------------------------------------------------------


            global initParams;

            if isfield(optimization_method, 'Particle_Swarm_Optimization') ...
                    && optimization_method.Particle_Swarm_Optimization

                [theta, objective] = OptAlgorithms.Particle_Swarm_Optimization(opt, params);

            elseif isfield(optimization_method, 'Differential_Evolution') ...
                    && optimization_method.Differential_Evolution

                [theta, objective] = OptAlgorithms.Differential_Evolution(opt, params);

            elseif isfield(optimization_method, 'Metropolis_Adjusted_Differential_Evolution') ...
                    && optimization_method.Metropolis_Adjusted_Differential_Evolution

                [theta, objective] = OptAlgorithms.Metropolis_Adjusted_Differential_Evolution(opt, params);

            elseif isfield(optimization_method, 'Markov_Chain_Monte_Carlo') ...
                    && optimization_method.Markov_Chain_Monte_Carlo

                [theta, objective] = OptAlgorithms.Markov_Chain_Monte_Carlo(opt, params);

            elseif isfield(optimization_method, 'Deterministic_algorithm_fmincon') ...
                    && optimization_method.Deterministic_algorithm_fmincon


                loBound = opt.paramBound(:,1);
                upBound = opt.paramBound(:,2);

                options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter',...
                    'TolX',1e-6,'TolCon',1e-6,'TolFun',1e-6,'MaxIter',500);

                try
                    [theta, objective, exitflag, output, ~, grad] = fmincon( OptAlgorithms.FUNC, ...
                        initParams, [],[],[],[], loBound, upBound, [], options);
                catch exception
                    disp('Errors in the MATLAB build-in optimizer: fmincon. \n Please check your input parameters and run again. \n');
                    disp('The message from fmincon: %s \n', exception.message);
                end

            else

                error('The method you selected is not provided in this programme \n');

            end

        end

        function mutant_struct = discreteMutation(opt, structure)
% -----------------------------------------------------------------------------
% This is the mutation part of the upper-level structure optimization algorithm
% 
% First of all, two random selected structures are prepared; 
% Then the optimal structure until now is recorded;
% Lastly, the mutant_struct = rand_struct_1 &+ rand_struct_2 &+ optima_struct
% -----------------------------------------------------------------------------


            % Record the optimal structure so far and select two random structures
            rand_struct_1 = structure(randi(structNumber), 1:opt.nZone);
            rand_struct_2 = structure(randi(structNumber), 1:opt.nZone);

            [~, id] = min(structure(:,opt.nZone+1));
            optima_struct = structure(id, 1:opt.nZone);

            % Preallocate of the mutation structure
            mutant_struct = zeros(1,nZone);

            for i = 1:opt.nZone
                if rand < 0.33
                    mutant_struct(i) = rand_struct_1(i);
                elseif (0.33 <= rand) && (rand<= 0.67)
                    mutant_struct(i) = rand_struct_2(i);
                else
                    mutant_struct(i) = optima_struct(i);
                end
            end

        end

        function trial_struct = discreteCrossover(opt, structure, mutant_struct)
% -----------------------------------------------------------------------------
% The crossover part of the upper-level structure optimization algorithm
%
% The generation of a new trial structure is achieved by randomly combining the original
% structure and mutation structure.
% -----------------------------------------------------------------------------


            trial_struct = zeros(1, opt.nZone);

            for i = 1:opt.nZone
                if rand < 0.5
                    trial_struct(i) = structure(i);
                else
                    trial_struct(i) = mutant_struct(i);
                end
            end

        end

    end % upper level


end
% =============================================================================
%              The MATLAB library for optimization case studies
% 
%      Copyright © 2015-2016: Qiaole He
% 
%      Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
% 
%  All rights reserved. This program and the accompanying materials
%  are made available under the terms of the GNU Public License v3.0 (or, at
%  your option, any later version) which accompanies this distribution, and
%  is available at http://www.gnu.org/licenses/gpl.html
% =============================================================================
