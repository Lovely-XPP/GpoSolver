function D = degree(P, vars)

if (nargin == 1)
    vars = symvar(P);
else
    vars = vars(:)';
end

s = size(P);
p = P(:);
dg = zeros(length(p),1);

for i = 1:length(p)
    dg(i) = deg(p(i), vars);
end

D = reshape(dg, s);

end

function dg = deg(p, vars)
 cv = char(vars);
 if (numel(vars) > 1)
    cv = cv(2:end-1);
 elseif (numel(vars) == 0)
     dg = 0;
     return;
 end
 
 cp = char(p);
 
 dg = feval(symengine, 'degree', cp);

if (numel(vars) == 0)
    dg = 0;
    return;
end

if (isempty(symvar(p)))
    dg = 0;
    return;
end

[~, mons] = coeffs(p,vars);

dg = zeros(size(mons));
for i = 1:numel(mons) 
    
    for j = 1:numel(vars)
        [~, m] = coeffs(mons(i), vars(j));
        dg(i) = dg(i) + length(sym2poly(m(1))) - 1;
    end
end
dg = max(dg);
end
