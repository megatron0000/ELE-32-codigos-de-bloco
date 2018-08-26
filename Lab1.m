origem = [0.5 0.2 0.1];
p = 0.5;
while  p(end)> 10^-6
   p = [p, p(end)/2]; 
end

chances = Test(p, @(x) Hamming47(x), @(x) DeHamming47(x));

loglog(1./p,chances, 1./p, p);