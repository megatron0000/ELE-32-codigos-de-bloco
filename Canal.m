function [ resultado ] = Canal( x, p )
resultado = mod(x + (rand(size (x,1),size(x,2)) < p), 2);
end

