function relax = gporelax(prob, order)
% GPORELAX - construct a LMI problem class relaxation 
%
% Given problem class parametrization PROB and relaxation order ORDER,
% constructs LMI relaxation of PROB or order ORDER.
%
% RELAX = GPORELAX(PROB, ORDER)
%
% If the parameter ORDER is not given, LMI relaxation of the lowest
% possible order will be constructed. The function returns the LMI
% relaxation RELAX.
%
%   See also GPOGENERATOR, GPOSOLVE.

% (c) J. Heller, 2014-2015

% Test the inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isfield(prob, 'cons'))
    prob.cons = [];
end

prob = default_param(prob, 'vars', []);
prob = default_param(prob, 'rpars', []);
prob = default_param(prob, 'ppars', []);
prob = default_param(prob, 'obj', []);
prob = default_param(prob, 'cons', []);

no_cs = numel(prob.cons);
no_vars = numel(prob.vars);

if (no_vars == 0)
    error('Problem does not contain any variables');
end

if (isempty(prob.cons) && isempty(prob.obj))
    error('Problem does not contain neither "obj" nor "cons"');
end

% Symbolic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prob.vars = prob.vars(:)';
cv = char(prob.vars);
if (no_vars > 1)
    cv = cv(9:(end-2));
end

fprintf(1, 'Optimization variables ... %s \n', cv);

% Is the cost function a polynomial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~ispoly(prob.obj, prob.vars))
    error('Cost function is not a polynomial');
end

% Convert constraints to the canonical form %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ccons = cell(0,0);
tvars_map = containers.Map('KeyType','char', 'ValueType','any');

if (~iscell(prob.cons))
    prob.cons = num2cell(prob.cons, 1);
end

for i = 1:no_cs
    con = prob.cons{i};
    no_cons = numel(con);
    siz = size(con);
    
    lhs = sym(zeros(no_cons, 1));
    rhs = sym(zeros(no_cons, 1));
    
    mop = '';
    for j=1:no_cons
        cstr = char(con(j));

        op = char(feval(symengine, 'op', cstr, 0));
        if (strcmp(mop, ''))
            mop = op;
        elseif (~strcmp(mop, op))
            error(['Contraints operators not consistent: ' cstr]);
        end
        
        lhs(j) = feval(symengine, 'op', cstr, 1);
        rhs(j) = feval(symengine, 'op', cstr, 2);
    end
    
    lhs = reshape(lhs, siz);
    rhs = reshape(rhs, siz);
        
    if (~ispoly([lhs, rhs]))
        error(['Not a polynomial constraint: ' cstr]);
    end
    
    if (max(siz) > 1 && min(siz) == 1)
        % linear constraint
        if (strcmp(op, '_leequal'))
            for c = 1:max(siz)
                [ccons{end + 1}, tvars_map] = poly2tpoly(rhs(c) - lhs(c), prob.vars, tvars_map); %#ok<AGROW>
            end
        elseif (strcmp(op, '_SYM_equal'))
            for c = 1:max(siz)
                [ccons{end + 1}, tvars_map] = poly2tpoly(lhs(c) - rhs(c), prob.vars, tvars_map); %#ok<AGROW>
                [ccons{end + 1}, tvars_map] = poly2tpoly(rhs(c) - lhs(c), prob.vars, tvars_map); %#ok<AGROW>
            end
        else
            error(['Binary operator not recognized: ' cstr '\n' ...
                'Only operators <=, >=, == are supported']);
        end
    elseif (siz(1) == siz(2))
        % scalar or LMI/PMI constraint
        if (strcmp(op, '_leequal'))
            [ccons{end + 1}, tvars_map] = poly2tpoly(rhs - lhs, prob.vars, tvars_map); %#ok<AGROW>
        elseif (strcmp(op, '_SYM_equal'))
            if (siz(1) == 1)
                % test constraint as a candidate for moment substitution
                % TODO
            end
            [ccons{end + 1}, tvars_map] = poly2tpoly(lhs - rhs, prob.vars, tvars_map); %#ok<AGROW>
            [ccons{end + 1}, tvars_map] = poly2tpoly(rhs - lhs, prob.vars, tvars_map); %#ok<AGROW>
        else
            error(['Binary operator not recognized: ' cstr '\n' ...
                   'Only operators <=, >=, == are supported']);
        end
    end
end

prob.ccons = ccons;
prob.tvars_map = tvars_map;

% Relaxation order %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, 'Computing problem polynomial degree ... ');

dc = [];
for i = 1:numel(prob.ccons)
    dc = [dc, max(degree(prob.ccons{i}, prob.vars))]; %#ok<AGROW>
end

do = degree(prob.obj, prob.vars);
prob.max_degree = max([do, dc]);
prob.rankshift = max(1, ceil(max(dc)/2));

fprintf(1, '%d\n', prob.max_degree);

% LMI problem construction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (prob.max_degree == 1)
    % LMI problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf(1, 'Building LMI problem ... ');
    no_cs = numel(prob.ccons);
    %prob.Cs = cell(no_cs, 1);

    blk_sz = cellfun(@(x) numel(x), prob.ccons, 'UniformOutput', true);
    [srt, idx] = sort(blk_sz);
    lp_pos = find(srt <= 1, 1, 'last');
    
    if ( lp_pos > 2)
        % we can construct an LP block
        prob.Cs = prob.ccons(idx((lp_pos + 1):end));
        prob.Cs{end + 1} = diag([prob.ccons{idx(1:lp_pos)}]);
       % prob.ccons = prob.Cs;
    else
        prob.Cs = prob.ccons;
    end
    
    prob.order = 1;
    
    fprintf(1, 'done\n');
else
    % LMI relaxation problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf(1, 'Computing minimum relaxation degree ... ');
    
    rc = ceil(dc/2);
    ro = ceil(do/2);
    prob.min_order = max([1, rc, ro]);
    fprintf(1, '%d\n', prob.min_order);
    
    if (nargin == 1)
        order = prob.min_order;
    elseif (prob.min_order > order)
        order = prob.min_order;
        warning(['Requested relaxation order is less that the minimal possible order. '...
            'Building relaxation of order ' sprintf('%d', order)]);
    end
    
    prob.order = order;
    
    % Moment matrix

    fprintf(1, 'Computing moment matrix ... ');
    
    no_vars = numel(prob.vars);
    no_cs = numel(prob.ccons);
    Cs = cell(no_cs + 1, 1);    
    
    b = basis(prob.vars, order);
    MM = b' * b;
    Cs{1} = MM;
    
    if (isempty(prob.obj))
        prob.obj = trace(MM);
    end
        
    fprintf(1, 'done\n');

    % Localizing matrices

    fprintf(1, 'Computing localizing matrices ... ');

    for i = 1:no_cs
        o = uint32(order - rc(i));
        sz = nchk(no_vars + o, o);
        Cs{i + 1} = expand(kron(MM(1:sz, 1:sz), prob.ccons{i}));
    end

    prob.Cs = Cs;
    fprintf(1, 'done\n');
end

% Identify monomials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, 'Identifying monomials ... ');

no_cons = numel(prob.Cs);
mons_idx = containers.Map();

if (prob.max_degree == 1)
    mons_sym = basis(prob.vars, 1);
else
    mons_sym = basis(prob.vars, prob.order * 2);
end

mons_str = arrayfun(@(x) char(x), mons_sym, 'UniformOutput', false);

no_mons = numel(mons_str);
prob.Fs = cell(no_cons, no_mons);

cs_sizes = zeros(no_cons, 1);
for i = 1:no_cons
    [cs_sizes(i), ~] = size(prob.Cs{i});
end

for j = 1:no_mons
    mons_idx(mons_str{j}) = j;
    for i = 1:no_cons
        prob.Fs{i,j} = sssmatrix(cs_sizes(i));
    end
end

fprintf(1, '%d\n', numel(mons_sym));

% LMI decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, 'Decomposing LMI constraints ... ');

for l = 1:no_cons
    for i = 1:cs_sizes(l)
        for j = i:cs_sizes(l)
            [c, m] = coeffs(prob.Cs{l}(i,j), prob.vars);
            for k = 1:numel(m)
                prob.Fs{l, mons_idx(char(m(k)))}.push(i, j, c(k));
            end            
        end
    end
end

fprintf(1, 'done\n');

% Cost function decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, 'Decomposing objective function ... ');

[c, m] = coeffs(prob.obj, prob.vars);
cobj = sym(zeros(no_mons, 1));
for k = 1:numel(m)
    cobj(mons_idx(char(m(k)))) = c(k);
end

% Cut off the first constant monomial 
prob.cobj_const = cobj(1);
prob.cobj = cobj(2:end);
prob.mons_sym = mons_sym(2:end);
prob.cs_sizes = cs_sizes;

% Rename prob to relax
relax = prob;

fprintf(1, 'done\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exp = default_param(exp, param_name, param_value)
    if (~isfield(exp, param_name))
        exp = setfield(exp, param_name, param_value);
    end
end

function [tp, tmp_map] = poly2tpoly(p, vars, tmp_map)
siz = size(p);
tp = sym(zeros(size(p,1), size(p,2)));

for j = 1:numel(p)
    no_tvars = tmp_map.Count;
    [c, t] = coeffs(p(j), vars);
    
    no_coeffs = numel(c);
    
    c_tmp = sym(zeros(no_coeffs, 1));
    
    for i = 1:no_coeffs
        if (isnumber(c(i)))
            c_tmp(i) = c(i);
        elseif (ispoly(c(i)))
            tvar_str = sprintf('aux%04d', no_tvars);
            no_tvars = no_tvars + 1;
            tvar = sym(tvar_str, 'real');
            
            c_tmp(i) = tvar;
            tmp_map(tvar_str) = c(i);
        else
            error(['Coefficient is not a polynomial: ' char(c(i))]);
        end
    end
    
    tp(j) = dot(c_tmp, t);
end
end

function B = basis(vars, deg)
no_vars = numel(vars);
P = pbasis(no_vars, deg);
no_mons = size(P, 1);
B = sym(zeros(1, no_mons));

for i = 1:no_mons
    B(i) = prod(vars.^P(i,:));
end
end



