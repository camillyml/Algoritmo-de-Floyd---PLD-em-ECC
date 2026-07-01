# Algoritmo-de-Floyd---PLD-em-ECC
# Colisão em Curvas Elípticas usando Floyd e BSGS

Implementação de uma prova de conceito para encontrar dois inteiros distintos `x` e `y` tais que `xG = yG` em curvas elípticas com `a = b = 2` sobre corpos finitos de 25, 30, 34 e 40 bits.

## Algoritmos
- **Método simples (tartaruga e lebre)**: tartaruga avança +1, lebre dobra.
- **Variação eficiente (baby-step giant-step)**: busca em `O(sqrt(p))`.
