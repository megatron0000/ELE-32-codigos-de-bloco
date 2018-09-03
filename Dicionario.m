function [ dict ] = Dicionario( Ht , pesoMaximo)
% DICIONARIO[Ht] gera um mapa de síndromes para erros,
% em que o erro é o de menor peso de Hamming possível.
% Mais especificamente, é de representações em string
% de síndromes para erros.
% `pesoMaximo` define o limite para o peso dos erros que
% serão incluídos no mapa

% Tamanhos da palavra-código e da síndrome
tamErro = size(Ht, 1);

dict = containers.Map();

% Variáveis de loop
tamSubconjunto = 1;
erro = zeros(1, tamErro);

while tamSubconjunto <= pesoMaximo
    triplasIndices = nchoosek(1:tamErro, tamSubconjunto);
    for i=1:nchoosek(tamErro, tamSubconjunto)
        sindrome = mod(sum(Ht(triplasIndices(i,:), :), 1) , 2);
        sindromeHash = num2str(sindrome);
        if dict.isKey(sindromeHash)
            continue
        end
        erro(triplasIndices(i,:)) = 1;
        dict(sindromeHash) = erro;
        erro(triplasIndices(i,:)) = 0;
    end
    tamSubconjunto = tamSubconjunto + 1;
end

end