#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <boost/dynamic_bitset.hpp>
#include <boost/container/vector.hpp>
#include <boost/unordered_map.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/niederreiter_base2.hpp>
#include <./Timer.cpp>

typedef boost::dynamic_bitset<> Vetor;
typedef boost::container::vector<Vetor> Matriz;

/*
*
* Ht: Matriz de verificação de paridade original
*
* Ht_long: Referência de onde guardar o array em que cada elemento é uma linha de Ht
* 
* Erros_1: Referência de onde guardar o
    array em que cada elemento é um erro de peso=1 (mas só a porção de informação desse erro)
*   exemplo: se há 5 bits na palavra-código e os 2 primeiros são de informação, então um elemento 
*   de Erros_1 pode ser 3=0b11, significando erro em cada um dos 2 bits de informação
*
* n_linhas, n_colunas, n_informacao: Quantas linhas e colunas tem Ht, e quantos bits são de informação
*   na palavra código (assume-se que são correspondem às primeiras linhas de Ht)
*/
void __gera_auxiliares(
    const Matriz& Ht,
    unsigned long* Ht_long,
    unsigned long* Erros_1,
    int n_linhas,
    int n_colunas,
    int n_informacao
) {
    /**
    * Só necessários os ulongs correspondentes a cada síndrome
    */
    for (int i = 0; i < n_linhas; ++i)
    {
        Ht_long[i] = Ht[i].to_ulong();
    }

    /** 
    * Porção de informação dos erros associados a cada síndrome
    * (por isso só importam os erros em bits de informação)
    */
    for (int i = 0; i < n_linhas; ++i)
    {
        // Erros depois de n_informacao linhas correspondem somente a erros de bit
        // de paridade, portanto aparecem como 0 já que só se coletam os erros de informação
        Erros_1[i] = i < n_informacao ? ((unsigned long)1) << (n_informacao - i - 1) : 0;
    }
}





/**
* input_csv: Espera-se formato de valores 0 ou 1 separados por vírgula e \n, sem espaços
*
* M: Matriz onde armazenar a leitura. 
*   Restrição: linha/coluna pode ter, no máximo (inclusive), 64 bits (cabe num unsigned long)
*/
void carrega_matriz(
    std::ifstream& input_csv,
    Matriz& M
) {
    std::string numero_string;

    while (getline(input_csv, numero_string)) {
    
        int commaPos;
        
        while ((commaPos = numero_string.find(",")) != std::string::npos)
            numero_string.erase(commaPos, 1);

        Vetor numero (numero_string.length(), stoul(numero_string, 0, 2));

        // std::cout << numero << "\n";

        M.push_back(numero);

    }
}

int testa_carrega_matriz() {
    std::string Htcsv = "./dados/Ht.csv";
    std::ifstream input {Htcsv};

    // if(!input)
    //     std::error("could not open " + Htcsv);

    std::string numero_string;
    Matriz Ht;

    while (getline(input, numero_string)) {
    
        int commaPos;
        
        while ((commaPos = numero_string.find(",")) != std::string::npos)
            numero_string.erase(commaPos, 1);

        Vetor numero (numero_string.length(), stoul(numero_string, 0, 2));

        std::cout << numero << "\n";

        Ht.push_back(numero);

    }

    std::cout << "Finished loading. Now replicating:" << std::endl;

    std::for_each(Ht.begin(), Ht.end(), [](Vetor line) {
        std::cout << line << std::endl;
    });
}

/**
* Auxiliar para `popular_dict`
*
*/
void __popular_dict(
    boost::unordered_map<unsigned long, unsigned long>& dict,
    unsigned long* Ht_long,
    unsigned long* Erros_1,
    int n_linhas,
    int n_colunas,
    int n_informacao,
    unsigned long sindrome_parcial,
    unsigned long erro_parcial,
    int linha_inicio,
    int niveis_restantes
) {
    if (niveis_restantes == 1) {
        unsigned long sindrome, erro;
        for (int i = linha_inicio; i < n_linhas; ++i)
        {
            sindrome = sindrome_parcial ^ Ht_long[i];
            erro = erro_parcial ^ Erros_1[i];
            if (dict.find(sindrome) == dict.end())
                dict[sindrome] = erro;
        }
        return;
    }

    for (int i = linha_inicio; i < n_linhas - (niveis_restantes - 1) /*iterar todas as síndromes*/; ++i)
    {
        __popular_dict(
            dict,
            Ht_long,
            Erros_1,
            n_linhas,
            n_colunas,
            n_informacao,
            sindrome_parcial ^ Ht_long[i],
            erro_parcial ^ Erros_1[i],
            i+1,
            niveis_restantes-1
        );
    }

}

/**
* 
* dict: Onde guardar o mapa "síndrome->erro"
*
* Ht: Matriz de verificação de paridade (transposta). Assume-se que suas últimas linhas
* correspondem aos bits de informação (identidade em cima, outras linhas embaixo)
*
* n_linhas, n_colunas, n_informacao: Quantas linhas e colunas tem Ht, e quantos bits são de informação
*   na palavra código (assume-se que são correspondem às primeiras linhas de Ht)
*
* peso_maximo: O peso dos maiores erros de informação a serem catalogados em `dict`
* 
*/
void popular_dict(
    boost::unordered_map<unsigned long, unsigned long>& dict,
    const Matriz& Ht,
    int n_linhas,
    int n_colunas,
    int n_informacao,
    int peso_maximo
) {

    unsigned long *Ht_long = new unsigned long[n_linhas];
    unsigned long *Erros_1 = new unsigned long[n_linhas];
    __gera_auxiliares(
        Ht,
        Ht_long,
        Erros_1,
        n_linhas,
        n_colunas,
        n_informacao
    );

    for (int i = 1; i <= peso_maximo; ++i)
    {
        __popular_dict(
        dict,
        Ht_long,
        Erros_1,
        n_linhas,
        n_colunas,
        n_informacao,
        0,
        0,
        0,
        i
    );
    }

}

void testa_popular_dict() {
    std::ifstream Htcsv("dados/Ht.csv");
    Matriz Ht;
    boost::unordered_map<unsigned long, unsigned long> dict;

    carrega_matriz(Htcsv, Ht);
    int n_linhas = Ht.size();
    int n_colunas =  Ht[0].size();
    int n_informacao = n_linhas - n_colunas;
    
    popular_dict(
        dict,
        Ht,
        n_linhas,
        n_colunas,
        n_informacao,
        3
    );

    assert(dict.find(0) == dict.end());
    assert(dict.at(1099511595008) == 34359738368); // 1ª linha de Ht -> erro só no 1º bit de info
    assert(dict.at(1064615018496) == 17179869184); // 2ª linha de Ht -> erro só no 2º bit de info
    assert(dict.at(34896576512) == 51539607552); // 1ª+2ª linhas de Ht -> erros nos 2 MSB de info
    assert(dict.at(1052568944640) == 55834574848); // 1ª+2ª+4ª linhas de Ht -> erros nos 1º,2º,4º bits de info
    assert(dict.at(549755813888) == 0); // erro composto só pela 37ª linha de Ht não aparece nos bits de info
    assert(dict.at(481038172167) == 1); // mas se erro for de 36ª+37ª linhas de Ht -> erro no último bit de info 
                                         // (e ignora o bit de paridade)
}


/**
* Retorna nova matriz em que cada linha é uma coluna da matriz original.
*/
Matriz& transposta(
    const Matriz& M
) {
    int antes_colunas = M[0].size();
    int antes_linhas = M.size();
    Matriz* resultado = new Matriz(antes_colunas);
    
    for (int i = 0; i < antes_colunas; ++i)
    {
        (*resultado)[i] = Vetor (antes_linhas);
        // Se não for Vetor&, boost copia implicitamente
        Vetor& linha = (*resultado)[i];
        for (int j = 0; j < antes_linhas; ++j)
        {
            linha[antes_linhas - j - 1] = M[j][antes_colunas - i - 1];
        };
    }
    return *resultado;
}

void testa_transposta() {
    Matriz m = Matriz(3);
    /* 
    * m == [1 1; 0 1; 0 0]
    */
    m[0] = Vetor(2, 3);
    m[1] = Vetor(2, 1);
    m[2] = Vetor(2, 0);

    Matriz mt = transposta(m);

    assert(mt[0] == Vetor(3, 4));
    assert(mt[1] == Vetor(3, 6));
}

/**
* Calcula M.v
*/
Vetor& mult(const Matriz& M, const Vetor& v) {
    int n_linhas = M.size(), n_colunas = v.size();
    Vetor* resultado = new Vetor(n_linhas);

    for (int i = 0; i < n_linhas; ++i)
    {
        (*resultado)[n_linhas - i - 1] =  (M[i] & v).count() & 1;
    }

    return *resultado;
}

/**
* Calcula M.N
*/
Matriz& mult(const Matriz& M, const Matriz& N) {
    int n_linhas = M.size();
    int n_colunas = N[0].size();
    Matriz pre_resultado = Matriz(n_colunas);

    Matriz Nt = transposta(N);

    for (int i = 0; i < n_colunas; ++i)
    {
        pre_resultado[i] = mult(M, Nt[i]);
    }

    Matriz& resultado = transposta(pre_resultado);

    return resultado;
}

void testa_mult() {
    // M == [1 0 1; 0 0 1]
    Matriz M = Matriz(2);
    M[0] = Vetor(3, 5);
    M[1] = Vetor(3, 1);

    // v == [0 1 1]
    Vetor v = Vetor(3, 3);

    // resultado deve ser [1 1]
    assert(mult(M, v) == Vetor(2, 3));

    // [1 1; 0 0; 1 1]
    Matriz N = Matriz(3);
    N[0] = Vetor(2, 3);
    N[1] = Vetor(2, 0);
    N[2] = Vetor(2, 3);

    assert(mult(M, N).size() == 2);
    assert(mult(M, N)[0] == Vetor(2, 0));
    assert(mult(M, N)[1] == Vetor(2, 7));
}

/**
* Modifica `Transmitido`, com chance `p` de inverter cada bit
* de cada palavra código.
*
* Transmitido: Assume-se que cada linha é uma palavra-código
*/
void canal(Matriz& Transmitido, double p) {
    int count=0;
    int n_linhas = Transmitido.size();
    int n_colunas = Transmitido[0].size();
    Vetor* linha;
    boost::random::niederreiter_base2 gen(4);
    boost::random::uniform_01<double> random;
    for (int i = 0; i < n_linhas; ++i)
    {
        linha = &Transmitido[i];
        for (int j = 0; j < n_colunas; ++j)
        {
            if (random(gen) < p) {
                (*linha)[j].flip();
                count++;
            }
        }
    }
    // std::cout << "p ≃ " << ((double) count) / (n_linhas*n_colunas) << std::endl;
}

/**
* amostras_informacao: Lista de palavras de informacao a serem enviadas
*   Espera-se uma palavra por linha, elementos separados por vírgulas (sem espaço)
* 
* Ht_csv: Matriz de verificação de paridade (transposta). Mesmo formato de `amostras_informacao_csv`
*
* Gt_csv: Matriz de geração do código (transposta). Mesmo formato de `amostras_informacao_csv`
*
* p: Lista de chances de o canal BSC inverter um bit transmitido
*
* peso_maximo_memorizado: Caso se encontre síndrome com peso maior que isto, ela não será corrigida
* 
* Retorna lista de chances de erro de bit (uma para cada valor de p)
*/
boost::container::vector<double> desempenho(
    std::ifstream& amostras_informacao_csv,
    std::ifstream& Ht_csv,
    std::ifstream& Gt_csv,
    const boost::container::vector<double>& p,
    int peso_maximo_memorizado
) {
    boost::container::vector<double>* resultado = new boost::container::vector<double>(p.size());

    Matriz Ht;
    carrega_matriz(Ht_csv, Ht);

    /**
    * Exemplo de chave e valor:
    * - chave: 0b1001 , valor: 0b101 significa uma síndrome [1, 0, 0, 1] com erro associado [1, 0, 1, ...]
    *   em que as reticências indicam a parte do erro concernente aos bits de paridade (não importam)
    */
    boost::unordered_map<unsigned long, unsigned long> dict;

    int n_linhas = Ht.size();
    int n_colunas = Ht[0].size();
    int n_informacao = n_linhas - n_colunas;


    popular_dict(
        dict,
        Ht,
        n_linhas,
        n_colunas,
        n_informacao,
        peso_maximo_memorizado
    );

    Matriz Info;
    carrega_matriz(amostras_informacao_csv, Info);

    int n_amostras = Info.size();

    Matriz Gt;
    carrega_matriz(Gt_csv, Gt);
    Matriz G = transposta(Gt);

    for (int i_p = 0; i_p < p.size(); ++i_p)
    {
        Matriz Transmitido = mult(Info, G);

        canal(Transmitido, p[i_p]);

        Matriz Sindromes = mult(Transmitido, Ht);

        Matriz Transmitido_informacao = Matriz(n_amostras);
        for (int i = 0; i < n_amostras; ++i)
        {
            Transmitido_informacao[i] = Vetor(n_informacao);
            for (int j = 0; j < n_informacao; ++j)
            {
                Transmitido_informacao[i][n_informacao - j - 1] = Transmitido[i][n_linhas - j - 1];
            }
            // assert(Transmitido_informacao.size() == Transmitido.size());
            // assert(Transmitido_informacao[0].size() == Info[0].size());
            // assert(Transmitido_informacao[i] == Info[i]);
        }

        int n_erros = 0;
        int incr=0;
        unsigned long sindrome;
        Vetor correcao_nula = Vetor(n_informacao, 0);
        Vetor correcao;
        for (int i = 0; i < n_amostras; ++i)
        {
            sindrome = Sindromes[i].to_ulong();
            if (dict.find(sindrome) == dict.end())
                correcao = correcao_nula;
            else 
                correcao = Vetor(n_informacao, dict.at(sindrome));
            incr = (Transmitido_informacao[i] ^ correcao ^ Info[i]).count();
            // if (incr != 0) {
            //     std::cout << correcao << ": " << Transmitido_informacao[i] << " versus " << Info[i] << std::endl;
            //     std::cout << "diferença " << (Transmitido_informacao[i] ^ Info[i]) << std::endl;
            //     std::cout << "Originais: " << Transmitido[i] << " versus " << mult(Gt, Info[i]) << std::endl;
            //     std::cout << "diferença: " << (Transmitido[i] ^ mult(Gt, Info[i])) << std::endl;
            //     std::cout << "síndrome: " << Sindromes[i] << std::endl;
            //     std::cout << std::endl;
            // }

            n_erros += incr;
        }

        (*resultado)[i_p] = n_erros / ((double) n_amostras * n_informacao);
        std::cout << "parcial(" << i_p << "): " << (*resultado)[i_p] << std::endl;
    }
    return *resultado;
}

int main(int argc, char** argv) {

    int arg_peso_maximo = std::stoi(argv[1]);
    std::string arg_resultados = argv[2];

    std::ifstream Htcsv ("dados/Ht.csv");
    std::ifstream amostrasInput("dados/amostra-informacao.csv");
    std::ifstream GtInput("dados/Gt.csv");
    std::ofstream Resultados(arg_resultados);
    std::ofstream P("dados/lista-de-p.csv");

    boost::container::vector<double> p = boost::container::vector<double>(0);
    double p0 = 0.5;
    while(p0 > 1 /((double) 1000000)) {
        p.push_back(p0);
        p0 *= 0.5;
    }

    boost::container::vector<double> des = desempenho(amostrasInput, Htcsv, GtInput, p, arg_peso_maximo);
    Resultados << des[0];
    P << p[0];
    for (int i = 1; i < des.size(); ++i)
    {
        Resultados << ",";
        Resultados << des[i];
        P << ",";
        P << p[i];
    }
    Resultados << std::endl;
    P << std::endl;

    return 0;
}
