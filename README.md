# Esqueminha-boltz/Boltz Co-folding Job Builder

Interface web (roda no seu navegador) para gerar os arquivos `.yaml` de
entrada do [Boltz](https://github.com/jwohlwend/boltz) (predição de
estrutura proteína-ligante / co-folding), sem precisar escrever YAML na mão.

## O que essa ferramenta gera (não roda o Boltz por si)

Ela só monta os arquivos de entrada certinhos e um `.zip` pronto para você
baixar e rodar. Três modos, em abas separadas:

1. **Proteína única + Ligante único** — um complexo.
2. **Screening** — 1 proteína contra uma lista de ligantes (ex: triagem
  virtual de uma biblioteca de moléculas contra um alvo).
3. **Multidocking** — várias proteínas x várias ligantes (todas as
  combinações).

Cada `.zip` gerado já vem com:

- os arquivos `.yaml` prontos para o `boltz predict`;
- um script (`run.sh` ou `run_all.sh`) com o comando certo para rodar tudo
  de uma vez;
- (nos modos 2 e 3) um `manifest.csv` relacionando cada job ao nome da
  proteína/ligante original, e um `compile_results.py` que junta os
  resultados (confiança + afinidade, se marcada) num único CSV já
  ranqueado — assim você não precisa abrir pasta por pasta.
- um `LEIA-ME.txt` com o passo a passo.

## Como instalar e abrir a interface

```bash
pip install -r requirements.txt
python app.py
```

Vai aparecer um endereço no terminal, algo como:

```
Running on local URL:  http://127.0.0.1:7860
```

Abra esse endereço no navegador. Preencha os campos da aba desejada e
clique em **Gerar arquivos**. Baixe o `.zip`.

## Como rodar o Boltz de fato (depois de baixar o .zip)

Isso precisa ser feito no computador/servidor onde o Boltz e a GPU estão
instalados (pode ser o mesmo computador, ou um cluster/servidor remoto —
nesse caso, envie o `.zip` para lá).

```bash
pip install boltz -U        # só na primeira vez
unzip boltz_....zip -d meu_job
cd meu_job
bash run.sh          # modo 1 (single)
# ou
bash run_all.sh       # modos 2 e 3 (screening / multidocking)
```

Os resultados ficam em `resultados/predictions/<nome_do_job>/`, com o
arquivo `.cif` da estrutura predita e um `.json` com os escores de
confiança (e de afinidade, se ativado).

Para screening/multidocking, depois que o Boltz terminar, rode:

```bash
python compile_results.py
```

Isso cria `resumo_resultados.csv`, já ordenado (por afinidade quando
disponível, senão por confidence_score) — a forma mais rápida de ver quais
ligantes/combinações tiveram o melhor resultado.

## Formato dos campos de entrada

- **Sequência da proteína**: sequência de aminoácidos em código de 1 letra
  (ex: `MKTAYIAKQ...`). Pode colar com quebras de linha, o app limpa
  automaticamente.
- **Ligante**: SMILES (mais comum) ou código CCD (ex: `ATP`, `SAH`) de um
  componente já conhecido do PDB.
- **Listas de proteínas/ligantes** (modos 2 e 3): uma entrada por linha, no
  formato `nome,valor` — por exemplo:
  
  ```
  aspirina,CC(=O)Oc1ccccc1C(=O)O
  ibuprofeno,CC(C)Cc1ccc(cc1)C(C)C(=O)O
  ```
  
  Ou, se preferir, envie um arquivo CSV com colunas `nome` e
  `smiles`/`sequencia`.

## Opções avançadas

- **Gerar MSA automaticamente (`--use_msa_server`)**: recomendado deixar
  ligado. Sem isso, o Boltz roda em "modo sequência única", o que reduz
  bastante a qualidade da predição, a não ser que você já tenha um MSA
  pré-computado.
- **Calcular afinidade de ligação**: ativa o módulo de afinidade do
  Boltz-2, que estima o quão forte é a ligação (útil para ranquear
  ligantes em um screening). Só funciona bem para ligantes pequenos
  (a documentação do Boltz recomenda até ~56 átomos pesados).

## Referência

Formato de entrada baseado na documentação oficial do Boltz:
https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md
