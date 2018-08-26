function [ resultado ] = Hamming47( x )
resultado = mod(x*[eye(4),[1 1 1; 1 0 1; 1 1 0; 0 1 1]], 2); 
end

