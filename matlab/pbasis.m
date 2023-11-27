function P = pbasis(no_vars, deg)
sz = nchk(no_vars + deg, deg);
P = zeros(sz, no_vars);

r = 1;
for i = 0:deg
    s = nchk(no_vars - 1 + i, i);
    P(r:(r + s - 1), :) = powers(no_vars, i);
    r = r + s;
end
end

function A = powers(no_vars, deg)
sz = nchk(no_vars - 1 + deg, no_vars - 1);
A = zeros(sz, no_vars);

gpow(no_vars, deg, 1, 1);

    function r = gpow(nv, dg, row, col)
        if (col == no_vars)
            A(row, col) = dg;
            r = 1;
        else
            if (dg > 0)
                r = 0;
                for i = dg:-1:0
                    rd = gpow(nv - 1, dg - i, row + r, col + 1);
                    A((row + r):(row + r + rd - 1), col) = i;
                    r = r + rd;
                end
            else
                r = 1;
            end
        end
    end
end