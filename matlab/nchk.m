function r = nchk(n, k)

if (k > n)
    r = 0;
    return;
end

if (k * 2 > n)
    k = n-k;
end

if (k == 0)
    r = 1;
    return;
end

r = n;
for i = 2:k
    r = r * (n-i+1);
    r = r / i;
end
