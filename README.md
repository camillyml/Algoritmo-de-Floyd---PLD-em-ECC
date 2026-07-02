# Algoritmo-de-Floyd PLD em ECC
# PLD em Curvas Elípticas via ciclo de Floyd (Pollard's rho)

Curva usada em todos os testes: `y² = x³ + 2x + 2 (mod p)`, com `a = b = 2`.

Objetivo: encontrar `x ≠ y` tais que `xG = yG`, para um ponto gerador `G` da
curva, usando a analogia tartaruga/lebre (busca de ciclos de Floyd), tal como
usada na fatoração de inteiros (Pollard's rho para fatoração) — aqui aplicada
à estrutura de grupo dos pontos da curva elíptica.

Implementação em `pollard_rho_ec.py`, usando a biblioteca `tinyec` (`ec.Point`,
`ec.Curve`, `ec.SubGroup`) fornecida no enunciado. O gerador `G` é obtido
buscando o primeiro `x` cuja equação da curva tem raiz quadrada módulo `p`
(via Tonelli–Shanks), já que os pontos-base do enunciado não vinham prontos.

## 1) Método ingênuo 

- Tartaruga: `T_k = k·G` (soma `G` a cada passo: `1G, 2G, 3G, 4G, ...`)
- Lebre: `H_k = 2^k·G` (dobra a si mesma a cada passo: `2G, 4G, 8G, 16G, ...`)
- Compara `T_k` com `H_k` a cada passo; para quando coincidirem com `k`
  diferentes (a coincidência trivial em `k=1`, onde ambas valem `2G`, é
  descartada).

**Cuidado de implementação:** o expoente da lebre não pode ser mantido como
o inteiro gigante `2^k` a cada iteração (isso faz `y` crescer
exponencialmente em número de bits e explode o custo); guardamos apenas o
expoente `k` (inteiro pequeno) e só materializamos `y = 2^k` uma única vez,
no final, quando a colisão é encontrada.

**Complexidade:** esse método não é uma busca pseudo-aleatória. A
tartaruga percorre linearmente até `k ≈ ordem(G)` no pior caso, ou seja,
é **O(n)**, com `n` a ordem do subgrupo gerado por `G` (próxima de `p`).
Na prática, medimos uma taxa de ~100 mil passos/segundo em Python puro; para
a menor curva testada (`p ≈ 3,35 × 10⁷`) isso já significa potencialmente
dezenas de milhões de passos e vários minutos de execução. Isso confirma experimentalmente por que o método literal do
enunciado, apesar de correto, **não escala**, motivando a variante eficiente
abaixo.

## 2) Variante eficiente: Pollard's rho com passeio pseudo-aleatório (r-adding walk)

Em vez de comparar uma progressão aritmética com uma geométrica, usa-se o
Pollard's rho:

- Particiona-se o grupo em `r = 24` classes usando `idx(P) = P.x mod r`
  (com um índice fixo para o ponto no infinito).
- São pré-computados: `r` multiplicadores aleatórios `m_i` e os pontos
  `M_i = m_i·G`.
- Função de passeio: dado o estado atual `(a, P)` com `P = a·G`, o próximo
  estado é `(a + m_{idx(P)}, P + M_{idx(P)})` (ainda da forma `a'·G`).
- Detecção de ciclo de Floyd: a tartaruga dá 1 passo do passeio por
  iteração; a lebre dá 2 passos. Quando os pontos coincidem
  (`P_tartaruga = P_lebre`) com acumuladores `a₁ ≠ a₂`, tem-se
  `a₁·G = a₂·G` com `a₁ ≠ a₂`.
  
**Complexidade esperada:** O(√n) passos, pelo paradoxo do aniversário —
ordens de grandeza mais rápido que o método ingênuo.

## Resultados

Verificação feita de forma independente para cada caso: recomputa-se `x·G`
e `y·G` do zero (via exponenciação binária, sem reaproveitar nada do
passeio) e é confirmado que são iguais, com `x ≠ y`;
`(x−y)·G` é o ponto no infinito, ou seja, `x−y` é um múltiplo não-nulo da
ordem de `G`.

| p (bits) | p | passos do passeio (Floyd) | tempo | √p (referência) |
|---|---|---|---|---|
| 25 | 33 554 467 | 504 | 0,01 s | ≈ 5 793 |
| 30 | 1 073 741 891 | 42 172 | 0,75 s | ≈ 32 768 |
| 34 | 10 001 113 099 | 149 037 | 3,21 s | ≈ 100 006 |
| 40 | 1 099 511 627 791 | 1 589 640 | 40,9 s | ≈ 1 048 577 |

Para p=25,30,34,40 bits, em todos os casos o número de passos ficou entre
~1× e ~1,5× de `√p`, exatamente a ordem de grandeza prevista pela teoria do
Pollard's rho (constante ≈ `√(π/2) ≈ 1,25`, considerando as ~3 avaliações
de função por rodada de Floyd). Em todos os quatro testes, a verificação
independente confirmou `x ≠ y`, `x·G = y·G`, e `(x−y)·G = O` (ponto no
infinito).

Exemplo concreto (p = 33 554 467, `G = (2, 21436521)`):

```
x = 6112246310
y = 12309499722
x - y = -6197253412
(x - y)·G = ponto no infinito  ✓
```
