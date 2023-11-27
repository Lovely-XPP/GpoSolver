function R = ispoly(P, vars)

if (nargin < 2)
    vars = symvar(P);
end

s = size(P);
p = P(:);
r = logical(zeros(length(p),1));

for i = 1:length(p)
    r(i) = is_poly(p(i), vars);
end

R = reshape(r, s);

end

function r = is_poly(p, vars)

r = true;

if (isempty(symvar(p)))
    return;
end

[c, t] = coeffs(p, vars);

if ((t(end) == 1) && ~isempty(intersect(symvar(c(end)), vars)))
    r = false;
end

end