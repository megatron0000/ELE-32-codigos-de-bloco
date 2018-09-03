function [ chancebit ] = Test( p, Codificacao, Decodificacao)
quantlinhas = 500000;
info = randi([0, 1], [quantlinhas, 4]);
chancebit = p;
codificado = Codificacao(info);
for i = 1:length(p)
    recebido = Canal(codificado, p(i));
    decodificado = Decodificacao(recebido);
    resultado = sum(sum(abs(decodificado - info)==1));
    chancebit(i) = resultado/(quantlinhas*3);
end
end