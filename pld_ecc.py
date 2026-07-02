# -*- coding: utf-8 -*-
"""
PLD em Curvas Elipticas

Curva: y^2 = x^3 + a*x + b (mod p), com a = b = 2.

Duas implementacoes:
  1. naive_tortoise_hare  -> versao do enunciado:
         tartaruga: 1G, 2G, 3G, 4G, ...   (T_k = k*G)
         lebre:     2G, 4G, 8G, 16G, ...  (H_k = 2^k*G)
     Serve de prova de conceito, mas e O(ordem de G) no pior caso.

  2. rho_efficient -> Pollard's rho: passeio pseudo-aleatorio
     (r-adding walk) sobre multiplos de G, com deteccao de ciclo de Floyd
     (tartaruga anda 1 passo do passeio, lebre anda 2 passos do mesmo
     passeio). Complexidade esperada O(sqrt(ordem de G)).
"""
import random
import time

import tinyec.ec as ec

def tonelli_shanks(n, p):
    """Retorna um y tal que y^2 = n (mod p), ou None se n nao for residuo
    quadratico mod p."""
    n %= p
    if n == 0:
        return 0
    if pow(n, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)

    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1

    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1

    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, t2i = 0, t
        while t2i != 1:
            t2i = (t2i * t2i) % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, (b * b) % p, (t * b * b) % p, (r * b) % p
    return r


def find_point(a, b, p, start_x=1):
    """Encontra o primeiro ponto (x, y) sobre y^2 = x^3+ax+b mod p, x>=start_x."""
    x = start_x
    while True:
        rhs = (x ** 3 + a * x + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None:
            return x, y
        x += 1


def build_curve(p, a=2, b=2, name=None):
    gx, gy = find_point(a, b, p)
    field = ec.SubGroup(p, (gx, gy), p, 1)
    curve = ec.Curve(a, b, field, name=name or f"E_{p}")
    assert not curve.is_singular(), "curva singular, escolha outros a,b"
    assert curve.on_curve(gx, gy)
    return curve, curve.g



# Metodo ingenuo: tartaruga soma G, lebre dobra

def naive_tortoise_hare(G, max_steps=None, report_every=2_000_000):
    T = G
    H = G
    x = 1        # T_k = x*G  (x cresce linearmente: 1,2,3,4,...)
    k = 1        # H_k = 2^k*G  (guardamos so o EXPOENTE k, nao 2^k!)
    steps = 0
    t0 = time.time()
    while True:
        steps += 1
        T = T + G
        x += 1
        H = H + H
        k += 1
        if T == H and x != k:
            y = 1 << k
            return {
                "x": x, "y": y, "steps": steps,
                "time": time.time() - t0,
            }
        if max_steps and steps >= max_steps:
            return None
        if report_every and steps % report_every == 0:
            print(f"    ... {steps} passos ({time.time()-t0:.1f}s), sem colisao ainda")

# Pollard's rho eficiente: r-adding walk + ciclo de Floyd
def rho_efficient(G, r=20, seed=None, max_steps=20_000_000):
    curve = G.curve
    p = curve.field.p
    rnd = random.Random(seed)

    # Tabela de "saltos" pseudo-aleatorios pre-computada: para cada uma das
    # r particoes do grupo, um multiplicador m_i e o ponto M_i = m_i*G.
    mults = [rnd.randint(2, p - 2) for _ in range(r)]
    jumps = [m * G for m in mults]

    def idx_of(P):
        if isinstance(P, ec.Inf):
            return 0
        return P.x % r

    def step(a, P):
        i = idx_of(P)
        return a + mults[i], P + jumps[i]

    def restart():
        a0 = rnd.randint(2, p - 2)
        return a0, a0 * G

    a1, R1 = restart()
    a2, R2 = a1, R1
    steps = 0
    t0 = time.time()
    while True:
        steps += 1
        a1, R1 = step(a1, R1)
        a2, R2 = step(a2, R2)
        a2, R2 = step(a2, R2)

        if R1 == R2:
            if a1 != a2:
                return {
                    "x": a1, "y": a2, "steps": steps,
                    "time": time.time() - t0,
                }
            # colisao trivial (a1 == a2) -> reinicia o passeio
            a1, R1 = restart()
            a2, R2 = a1, R1

        if steps >= max_steps:
            return None


def verify(G, x, y):
    curve = G.curve
    xg = x * G
    yg = y * G
    ok = (xg == yg) and (x != y)
    d = x - y
    dG = d * G
    is_inf = isinstance(dG, ec.Inf)
    return ok, d, is_inf


if __name__ == "__main__":
    primes = [
        ("25 bits", 33_554_467),
        ("30 bits", 1_073_741_891),
        ("34 bits", 10_001_113_099),
        ("40 bits", 1_099_511_627_791),
    ]

    print("=" * 70)
    print("PARTE 1: metodo ingenuo (tartaruga soma G, lebre dobra) -- so na")
    print("curva menor, pois e O(ordem de G) no pior caso (~10^5 pontos/s).")
    print("Limite de passos abaixo para nao rodar indefinidamente; aumente")
    print("se quiser deixar rodar ate convergir de fato (pode levar minutos).")
    print("=" * 70)
    label, p = primes[0]
    curve, G = build_curve(p)
    print(f"\n[{label}] p={p}, curva: {curve}")
    print(f"G = {G}")
    res = naive_tortoise_hare(G, max_steps=5_000_000, report_every=1_000_000)
    if res is None:
        print("  Nao convergiu dentro do limite de passos (esperado: e O(n)).")
    else:
        ok, d, is_inf = verify(G, res["x"], res["y"])
        print(f"  Colisao encontrada em {res['steps']} passos ({res['time']:.2f}s)")
        print(f"  x={res['x']}  y=2^{res['y'].bit_length()-1}  (x != y)")
        print(f"  Verificacao independente x*G == y*G: {ok}")
        print(f"  (x-y)*G eh o ponto no infinito? {is_inf}")

    print()
    print("=" * 70)
    print("PARTE 2: Pollard's rho eficiente (r-adding walk, ciclo de Floyd)")
    print("=" * 70)
    for label, p in primes:
        curve, G = build_curve(p)
        print(f"\n[{label}] p={p}")
        print(f"  curva: {curve}")
        print(f"  G = {G}")
        res = rho_efficient(G, r=24, seed=1234, max_steps=20_000_000)
        if res is None:
            print("  Nao convergiu dentro do limite de passos.")
            continue
        ok, d, is_inf = verify(G, res["x"], res["y"])
        print(f"  Colisao encontrada em {res['steps']} passos do passeio "
              f"({res['time']:.2f}s)")
        print(f"  x={res['x']}")
        print(f"  y={res['y']}")
        print(f"  x != y: {res['x'] != res['y']}")
        print(f"  Verificacao independente x*G == y*G: {ok}")
        print(f"  x-y = {d}")
        print(f"  (x-y)*G eh o ponto no infinito (O)? {is_inf}  "
              f"-> (x-y) e MULTIPLO da ordem de G")
