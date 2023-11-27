function [R, Z] = isnumber(P)

if (nargout < 2)
    s = size(P);
    p = P(:);
    r = arrayfun(@(x) is_number(x), p);
    R = reshape(r, s);
else
    s = size(P);
    p = P(:);
    r = arrayfun(@(x) is_number0(x), p, 'UniformOutput', false);
    r = [r{:}];
    R = reshape(r(1,:), s); 
    Z = reshape(r(2,:), s);
end

end

function r = is_number(p)
    r = true;

    try 
        dummy = double(p); %#ok<NASGU>
    catch
        r = false;
    end
end

function r = is_number0(p)
    r = [true; false];

    try 
        dummy = double(p); %#ok<NASGU>
        if (dummy == 0);
            r(2) = true;
        end
    catch
        r(1) = false;
    end
end