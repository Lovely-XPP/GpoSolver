function [status, sols] = gposolve(relax, rvals, pvals, params)
% GPOSOLVE - solve POP/PMI problem instance using SeDuMi SDP solver
%
% Given a LMI relaxation RELAX and parameter values RVALS and PVALS,
% solves the original POP problem instance. If possible, GPOSOLVE will
% certify global optimality of the solution. GPOSOLVE needs SeDuMi toolbox
% as a prerequisite.
%
% [STATUS, SOLS] = GPOSOLVE(RELAX, RVALS, PVALS, PARAMS)
%
% Function GPOSOLVE takes four parameters. RELAX is a LMI relaxation of a
% POP/PMI problem class. RVALS is a l x m matrix of residual parameter values,
% such that 'l' is the number of residual blocks and 'm' is the number of
% residual parameters. The value order in the rows of RVALS matrix must correspond
% to the order of variables in the vector RPARS from the problem class
% definition. PVAL is a vector of problem parameter values, the value order
% must correspond to the order of vector PPARS from the problem class definition. 
% Data structure PARAMS can have the following optional parameters:
%
% verbose  : set to  1 to output SDP solver progress information and
%            solution information. Default value is 0;
% sdptol   : SDP feasibility tolerance. Default value is 1e-03.
% maxnorm  : Maximum value of the l2 -norm of the solution of Pδ . If the
%            l2 -norm of the solution is too large, it may indicate that
%            additional bound constraints should be enforced or the
%            problem/residual parameters scaled. Default value is 1e+06.
% restol   : Numerical tolerance of the solution feasibility test. 
%            Default value is 1e-03.
% ranktol  : Numerical tolerance of the matrix rank computation. 
%            Default value is 1e-03.
% pivtol   : Numerical tolerance of the Gaussian elimination pivoting of the
%            solution extraction algorithm. If set to a negative value, dynamic
%            tolerance selection algorithm is used instead. Default value is −1.
%
% If GPOSOLVE recovers any solutions, it will return a n x m matrix SOLS
% which rows correspond to n solution vectors. The value order of the solution
% vectors corresponds to the variable order in vector VARS from the
% problem class definition. Parameter STATUS can take one of the following
% values:
% 
% 0        : A solution set was successfully extracted and the
%            global optimality was certified.
% 1        : A solution set was successfully extracted, however,
%            global optimality could not be certified.
% 2        : LMI relaxation RELAX could be solved, but the
%            l2 -norm of the solution was too large. This may
%            indicate that additional bound constraints should be enforced.
% 3        : LMI relaxation RELAX could be solved, but the
%            solution extraction algorithm could not extract any solutions.
% 4        : LMI relaxation RELAX is infeasible.
%
%   See also GPOGENERATOR, GPORELAX.

% (c) J. Heller, 2014-2015
     
if (nargin < 4)
    params = [];
end

params = default_param(params, 'method', 'gposolver');
params = default_param(params, 'verbose', 0);

params = default_param(params, 'sdptol',  1e-03);
params = default_param(params, 'maxnorm', 1e+06);
params = default_param(params, 'restol',  1e-03);
params = default_param(params, 'ranktol', 1e-03);
params = default_param(params, 'pivtol', -1);

if (strcmp(params.method, 'gposolver'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GpoSolver + SeDuMi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute cost function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ((numel(relax.rpars) ~= size(rvals, 2)))
        error('Size of parameter "rvals" must be nxm, where m==numel(relax.rpars)');
    end
    
    if (numel(relax.ppars) ~= numel(pvals))
        error('The number of elements of "pvals" must be equal to the number of elements of relax.ppars');
    end

    relax.rvals = rvals;
    relax.pvals = pvals(:)';
    no_ovals = size(rvals, 1);

    if (no_ovals > 0)
        % Cost function is parametrized
        cost = zeros(size(relax.cobj));
        cost_const = 0;
        
        for i = 1:no_ovals
            cost = cost + double(subs(relax.cobj, [relax.rpars, relax.ppars], [rvals(i, :), pvals]));
            cost_const = cost_const + double(subs(relax.cobj_const, [relax.rpars, relax.ppars], [rvals(i, :), pvals]));
        end
    else
        % cost function is not parametrized
        cost = double(relax.cobj);
        cost_const = double(relax.cobj_const);
    end
        
    % Substitute to the temporary variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    no_pvals = numel(pvals);
    relax.tvals_map = containers.Map('KeyType','char', 'ValueType','double');

    if (numel(relax.ppars) ~= no_pvals)
        error('Size of parameter "cvals" must be equal to the size of relax.ppars');
    end
            
    if (no_pvals > 0)
        no_tvars = relax.tvars_map.Count; 
        tvars_str = keys(relax.tvars_map);
        tvars_exp = values(relax.tvars_map);
        
        for i = 1:no_tvars       
            relax.tvals_map(tvars_str{i}) = double(subs(tvars_exp{i}, relax.ppars, pvals));
        end
    end
        
    % Construct SDP in SeDuMi format %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    if (~isfield(relax, 'Fs'))
        error('Problem is missing LMI constraints');
    end
    
    [no_cons, no_fs] = size(relax.Fs);
    
    relax.sdp_b = -cost;
    
    relax.sdp_c = [];
    
    for i = 1:no_cons
        F = relax.Fs{i, 1}.get_real_matrix(relax.tvals_map);
        relax.sdp_c = [relax.sdp_c; F(:)];
    end    
    
    relax.sdp_A = [];
    for i = 1:no_cons
        A = sparse(relax.cs_sizes(i)^2, no_fs - 1);
        for j = 1:(no_fs - 1)
            F = relax.Fs{i, j + 1}.get_real_matrix(relax.tvals_map);
            A(:, j) = -F(:);
        end
        relax.sdp_A = [relax.sdp_A; A];
    end
    
    K.s = relax.cs_sizes;
        
    % Solve SDP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    status = 0; % success
    sols = [];
    
    sedumi_params.fid = 0; 
    if (params.verbose)
        sedumi_params.fid = 1;
    end
    
    if (isfield(params, 'sedumi_eps'))
        sedumi_params.eps = params.sedumi_eps;
    else
        sedumi_params.eps = 1-20;
    end        
    
    [relax.sdp_x, relax.sdp_y, sedumi_info] = sedumi(relax.sdp_A, relax.sdp_b, relax.sdp_c, K, sedumi_params);
    
    if (params.verbose)
        fprintf('\n');
    end
        
    relax.sdp_obj = - (relax.sdp_b' * relax.sdp_y - cost_const);
        
    % Test SeDuMi status %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (sedumi_info.numerr == 2)
        status = 4;
        if (params.verbose)
            warning('Failure when solving SDP problem');
        end
        return;
    end
    
    if (sedumi_info.pinf == 1)
        if (params.verbose)
            warning('Solvers suspects primal infeasibility');
            warning('SDP problem may be undounded');
            warning('Try to enforce additional bound constraints');            
        end
    end
    
    if (sedumi_info.dinf == 1)
        if (params.verbose)
            warning('Solver suspects dual infeasibility');
            warning('SDP problem may be infeasible');
            warning('Try to scale problem variables');
        end
    end
    
    % Check feasibility %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    relax.sdp_res = sdp_residual(relax);
    
    if (params.verbose)
        fprintf('Checking SDP feasibility...\n');
    end
            
    if (relax.sdp_res < -params.sdptol)
        status = 4;
        if (params.verbose)
            warning(['Infeasible SDP: sdp_res = ' sprintf('%g', relax.sdp_res) ... 
                ' ( < ' sprintf('%g', -params.sdptol) ')' ]);
        end
        return;
    elseif (relax.sdp_res < 0 && params.verbose)
        fprintf('  Marginally feasible SDP: residual = %g\n', relax.sdp_res);
    elseif (params.verbose)
        fprintf('  Feasible SDP\n');
    end
        
    % Check solution norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (params.verbose)
        fprintf(1, 'Checking norm of the SDP solution ...\n');
    end
    
    if (norm(relax.sdp_y) > params.maxnorm)
        status = 2;
        if (params.verbose)
            warning(['SDP problem may be unbounded (norm of solution > ' ...
                sprintf('%g', -params.maxnorm) ')']);
        end
        return;
    end
    
    % Check global optimality for one solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sol = relax.sdp_y(1:numel(relax.vars))';
    
    feas = test_poly(sol, relax, params);
    
    if (feas)
        % One globally optimal solution found
        if (params.verbose)
            fprintf(1, ' Found one globally optimal solution\n');
        end
        
        status = 0;        
        sols = sol;
        return;
    end
    
    % Check global optimality for multiple solutions %%%%%%%%%%%%%%%%%%%%%%
    
    if (params.verbose)
        fprintf('Checking ranks of moment matrices...\n');
    end
    
    z = relax.sdp_c - relax.sdp_A * relax.sdp_y;
    no_vars = numel(relax.vars);
    oldrank = 1;
    rankdiff = 0;
    
    gopt_flag = 0;

    mm_size = relax.cs_sizes(1);
    MM = reshape(z(1:(mm_size*mm_size)), mm_size, mm_size);
    
    for i = 1:relax.order
        sz = nchk(no_vars + i, no_vars);
        [r, U, S] = mrank(MM(1:sz, 1:sz), params.ranktol);

        if (params.verbose)
            fprintf('  Moment matrix of order %2d has size %2d and rank %2d \n', ...
                i, sz, r);
        end
        
        if (oldrank >= r)
            rankdiff = rankdiff + 1;
        else
            rankdiff = 0;
        end
        
        if (rankdiff >= relax.rankshift)
            gopt_flag = 1;
            break;
        end      
        
        oldrank = r;
    end
    
    if (~gopt_flag)
      status = 1;
    end
    
    if (params.verbose)
        if (gopt_flag)
            fprintf('  Rank conditions may ensure global optimality\n');
        else
            fprintf(' Rank conditions cannot ensure global optimality\n');
        end
    end
    
    % Extract solutions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if (params.verbose)
        fprintf('Extracting solutions...\n');
    end    
    
    % Cholesky factor of the moment matrix
    [r, U, S] = mrank(MM(1:sz, 1:sz), params.ranktol);
    U = U(:, 1:r) * diag(sqrt(S(1:r)));
    
    if (params.pivtol > 0)
        [W, w] = rref(U', params.pivtol);
        W = W';
    else
        cU = cond(U);
        tol = 1e-10;
        [W, w] = rref(U', tol);
        W = W';
        while ((tol < 1) && ((cond(W) / cU) >= 1e4))
            tol = tol * 5;
            [W, w] = rref(U', tol);
            W = W';
        end
    end
    
    % Recover multiplication matrices    
    N = cell(1, no_vars);
    i = 0; 
    fail_flag = 0;
    pb = pbasis(no_vars, relax.order * 2);
    bc = bincoeffs(no_vars, relax.order * 2);
    
    while (~fail_flag && (i < no_vars))
        i = i + 1;
        mpb = pb(w,:); 
        mpb(:,i) = mpb(:,i) + 1;
 
        mind = ones(size(w));
        for k = 1:no_vars
            mind = mind + bc(no_vars + 1 - k, 1 + sum(mpb(:, k:no_vars), 2));
        end
        
        fail_flag = any(mind > size(U,1));
        if (~fail_flag)
            N{i} = W(mind,:);
        end
    end
    
    if (fail_flag)
        status = 3;
        relax.sols = {};
        if (params.verbose)
            fprintf('  Incomplete basis - no solutions extracted');
            fprintf('  Global optimality cannot be ensured');
        end
        return;
    end
           
    % Compute common eigenvalues of multiplication matrices
    cfs = rand(no_vars, 1); 
    cfs = cfs / sum(cfs);
    no_sols = size(N{1}, 1);
    M = zeros(no_sols);
    %cfs = [0.3 0.7];
    for i = 1:no_vars
        M = M + cfs(i) * N{i};
    end
    
    % Schur decomposition of M
    [Q,T] = schur(M);
    if (isempty(T))
        status = 3;
        relax.sols = {};
        if (params.verbose)
            fprintf('  Incomplete basis - no solutions extracted');
            fprintf('  Global optimality cannot be ensured');
        end
        return;
    end
         
    % Retrieve optimal solutions
    tsols = cell(1, no_sols);
    for i = 1:no_sols
        tsols{i} = zeros(no_vars, 1);
        for j = 1:no_vars
            tsols{i}(j) = Q(:, i)' * N{j} * Q(:, i);
        end
    end
    
    % Test validity of retrieved solutions
    sols = zeros(0, no_vars);
    for i = 1:no_sols
        feas = test_poly(tsols{i}', relax, params);
        if (feas)
            sols(end + 1, :) = tsols{i}'; %#ok<AGROW>
        end
    end
    
    no_sols = size(sols, 1);
    if (no_sols == 0)
        sols = [];
        if (params.verbose)
            fprintf('No solutions could be extracted \n');
        end        
    else
        if (params.verbose)
            fprintf('Extracted %d global optimal solutions\n', no_sols);
            if (gopt_flag)
                fprintf('Global optimality certified numerically\n');
            end
        end    
    end
        
elseif (strcmp(params.method, 'gloptipoly'))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GloptiPoly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mset clear;
    
    % Define variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:numel(relax.vars)
        cmd = ['mpol  ' char(relax.vars(i)) ';'];
        eval(cmd);
        if (params.verbose == 1)
            disp(cmd);
        end
    end  
    
    % Define constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C = [];
    no_pvals = numel(pvals);
    
    if (numel(relax.ppars) ~= no_pvals)
        error('Size of parameter "cvals" must be equal to the size of relax.ppars');
    end    
    
    for i = 1:numel(relax.cons)
        if (numel(relax.cons{i}) > 1)
            error('GloptiPoly does not support matrix/vector constraints');
        end
        
        if (no_pvals > 0)
            Cs_str = char(subs(relax.cons{i}, relax.ppars, pvals));
        else
            Cs_str = char(relax.cons{i});
        end
        
        cmd = ['C = [C, ' Cs_str '];' ];
        eval(cmd);
        if (params.verbose == 1)
            disp(cmd);
        end
    end
        
    % Compute cost function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    obj_sym = 0;
    
    if (numel(relax.rpars) > 0)
        for i = 1:size(rvals, 1)
            obj_sym = obj_sym + subs(relax.obj, [relax.rpars, relax.ppars], [rvals(i, :), pvals]);
        end
    else
        obj_sym = relax.obj;
    end
    
    cmd = ['F = ' char(obj_sym) ';'];
    eval(cmd);
    if (params.verbose == 1)
        disp(cmd);
    end
    
    % Solve relaxation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pars.fid = 0;
    if (params.verbose)
       pars.fid = 1; 
    end
    mset(pars);
       
    if (~isfield(relax, 'order'))
        if (isempty(C))
            cmd = 'P = msdp(min(F));';
        else
            cmd = 'P = msdp(min(F), C);';
        end
    else
        if (isempty(C))
            cmd = ['P = msdp(min(F), ' sprintf('%d', relax.order) ');'];
        else        
            cmd = ['P = msdp(min(F), C,' sprintf('%d', relax.order) ');'];
        end
    end
    
    eval(cmd);

    if (params.verbose == 1)
        disp(cmd);
    end
    
    [status, relax.sdp_obj] = msol(P);
    
    if (status == 1)
        sols = [];
        no_vars = numel(relax.vars);
        for i = 1:no_vars
            cmd = ['sols = [sols, squeeze(double(' char(relax.vars(i)) '))];'];
            eval(cmd);
        end
    end
else
    error(['Unrecognized method type: ' params.method]);
end
    
end % solve_relax

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exp = default_param(exp, param_name, param_value)
    if (~isfield(exp, param_name))
        exp = setfield(exp, param_name, param_value);
    end
end

function res = sdp_residual(relax)
    % Compute dual residual of an SDP problem
    z = relax.sdp_c - relax.sdp_A * relax.sdp_y;
    res = inf;    
    offset = 0;
    
    for i = 1:numel(relax.cs_sizes)
       s2 = relax.cs_sizes(i) * relax.cs_sizes(i); 
       Z = reshape(z(offset + (1:s2)), relax.cs_sizes(i), relax.cs_sizes(i));
       res = min(res, min(eig(Z)));
       offset = offset + s2;
    end
end

function feasible = test_poly(sol, relax, params)
    feasible = 1;

    if (params.verbose)
        fprintf(1, 'Testing solution ...\n');
    end

    % Compute objective function value
    no_ovals = size(relax.rvals, 1);
    objc = subs(relax.obj, relax.vars, sol);

    obj_val = 0;
    if (no_ovals > 0)
        for j = 1:no_ovals
            obj_val = obj_val + subs(objc, [relax.rpars, relax.ppars], [relax.rvals(j, :), relax.pvals]);
        end
    else
        obj_val = objc;
    end

    if (abs(obj_val - relax.sdp_obj) < params.restol)
        if (params.verbose)
            fprintf(1, '  Solution reaches SDP objective\n');
        end
    else
        if (params.verbose)        
            fprintf(1, '  Solution does not reach SDP objective\n');
        end
        feasible = 0;
       % return;
    end
       
    % Compute constraints values
    for i = 1:numel(relax.ccons)
        con = relax.ccons{i};
        pars = setdiff(symvar(con), relax.vars);

        if (~isempty(pars))
            pars_vals = arrayfun(@(x) relax.tvals_map(char(x)), pars);
            con = subs(con, pars, pars_vals);
        end
        
        con_val = double(subs(con, relax.vars, sol));
        if (numel(con_val) > 1)
            con_val = min(eig(con_val));
        end;
                
        if (con_val < -params.restol)
            if (params.verbose)
                fprintf(1, '  Infeasible constraint found\n');
            end
            feasible = 0;
            return;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bin = bincoeffs(no_vars, deg)
    bin = zeros(no_vars,deg + 2);
    bin(1,:) = 0:(deg + 1);

    for k = 2:no_vars
       bin(k, :) = cumsum(bin(k - 1, :));
    end
end

function [r, U, S] = mrank(M, tol)
    if (nargin < 2)
        tol = 1e-03;
    end
    
    [U, S] = svd(M); 
    % Detect rank drop
    S = diag(S); 
    n = length(S);
    
    drop = find(S(2:n) ./ S(1:n-1) < tol);
    if ~isempty(drop)
        r = drop(1);
    else
        r = n;
    end
end