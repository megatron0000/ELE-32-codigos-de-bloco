function [ resultado ] = DeHamming47( x )
resultado = zeros(size(x,1), 4);
Ht = [[1 1 1; 1 0 1; 1 1 0; 0 1 1];eye(3)];
sindrome = mod(x * Ht, 2);
for linha = 1:size(sindrome,1)
    aux = repmat(sindrome(linha, :), [7, 1]);
    aux = mod(aux - Ht, 2);
    aux = sum(aux, 2) == 0;
    resultado(linha, :) = mod(x(linha,1:4) + aux(1:4)', 2);
end
end