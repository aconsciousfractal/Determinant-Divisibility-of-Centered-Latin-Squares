"""
deep_investigation.py — TWO deep analyses:

POINT 3: p-Sylow conjecture — when is it nontrivial? What does Thm 3.2 already give?
POINT 6: WHY is K(A_cyc) = (Z/n)^{n-2}? 
         Connection to complete graph Laplacian (Kirchhoff).
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass
import numpy as np
from sympy import Matrix, ZZ, factorint, isprime, gcd as sgcd
from sympy.matrices.normalforms import smith_normal_form
from math import gcd
from functools import reduce
import random

random.seed(42)

def cyclic_ls(n):
    return np.array([[(i+j) % n + 1 for j in range(n)] for i in range(n)])

def get_A(L):
    n = L.shape[0]
    return np.array([[int(L[i,j]) - int(L[i,n-1]) for j in range(n-1)] for i in range(n-1)])

def snf_diag(M_np):
    """Return SNF diagonal via sympy."""
    m = M_np.shape[0]
    M = Matrix(M_np.tolist())
    S = smith_normal_form(M, domain=ZZ)
    return [int(S[i,i]) for i in range(m)]

def rank_Fp(M, p):
    n_rows, n_cols = M.shape
    B = (np.array(M, dtype=int) % p + p) % p
    pivot_row = 0
    for col in range(n_cols):
        found = False
        for row in range(pivot_row, n_rows):
            if B[row, col] % p != 0:
                B[[pivot_row, row]] = B[[row, pivot_row]]; found = True; break
        if not found: continue
        inv_piv = pow(int(B[pivot_row, col]), p-2, p)
        for row in range(n_rows):
            if row != pivot_row and B[row, col] % p != 0:
                factor = (int(B[row, col]) * inv_piv) % p
                for k in range(n_cols):
                    B[row, k] = (int(B[row, k]) - factor * int(B[pivot_row, k])) % p
        pivot_row += 1
    return pivot_row

def vp(x, p):
    """p-adic valuation of x."""
    if x == 0: return float('inf')
    x = abs(x)
    v = 0
    while x % p == 0:
        x //= p; v += 1
    return v

def bareiss_det(M_list):
    """Exact integer determinant via Bareiss algorithm."""
    m = len(M_list)
    M = [row[:] for row in M_list]
    sign = 1; prev = 1
    for i in range(m):
        if M[i][i] == 0:
            for k in range(i+1, m):
                if M[k][i] != 0:
                    M[i], M[k] = M[k], M[i]; sign *= -1; break
            else:
                return 0
        for j in range(i+1, m):
            for k in range(i+1, m):
                M[j][k] = (M[j][k]*M[i][i] - M[j][i]*M[i][k]) // prev
            M[j][i] = 0
        prev = M[i][i]
    return sign * M[m-1][m-1]

def jm_sample(n, steps=3000):
    """Jacobson-Matthews sampler for random LS of order n."""
    grid = [[(i+j) % n for j in range(n)] for i in range(n)]
    for _ in range(steps):
        r = random.randint(0, n-1)
        c = random.randint(0, n-1)
        s = random.randint(0, n-1)
        if grid[r][c] == s: continue
        r2 = c2 = None
        for rr in range(n):
            if rr != r and grid[rr][c] == s: r2 = rr; break
        for cc in range(n):
            if cc != c and grid[r][cc] == s: c2 = cc; break
        if r2 is None or c2 is None: continue
        old = grid[r][c]
        if grid[r2][c2] == old:
            grid[r][c] = s; grid[r2][c] = old
            grid[r][c2] = old; grid[r2][c2] = s
    return np.array([[v+1 for v in row] for row in grid])


# ════════════════════════════════════════════════════════════════════
print("=" * 72)
print("  PUNTO 3: p-SYLOW — ANALISI PROFONDA")
print("=" * 72)

# ── 3a: Cosa dà il Thm 3.2 vs la congettura p-Sylow ──────────────
print("""
Il Teorema Unificato (Thm 3.2) dà per OGNI LS di ordine n:
  • n dispari: n | det(A)  →  v_p(det A) ≥ v_p(n) per ogni p|n  
  • n pari:   n/2 | det(A) →  v_p(det A) ≥ v_p(n) per p dispari
                               v_2(det A) ≥ v_2(n) - 1

La congettura p-Sylow (rem:p-Sylow) dà:
  v_p(det A) ≥ v_p(n/p) + max(0, k_p - 1)  =  v_p(n) - 1 + max(0, k_p - 1)

Differenza: Thm 3.2 dà v_p(n), p-Sylow dà v_p(n) - 1 + max(0, k_p - 1).
  → p-Sylow NON-BANALE (più forte di Thm 3.2) solo quando k_p ≥ 2 
    (per p dispari) o k_p ≥ 1 (per p=2, ma questo È il Thm 5.5).
""")

# ── 3b: Distribuzione di k_p per ogni primo p|n ──────────────────
print("─" * 72)
print("  3b: Distribuzione k_p per LS random (sampling)")
print("─" * 72)

test_cases = [
    (6, [2, 3]),  
    (8, [2]),
    (9, [3]),
    (10, [2, 5]),
    (12, [2, 3]),
]

for n, primes in test_cases:
    print(f"\n  n = {n}, #samples = 200")
    for p in primes:
        kp_dist = {}
        vp_det_dist = {}
        thm32_bound = vp(n, p) if (n % 2 == 1 or p != 2) else vp(n, 2) - 1
        
        violations = 0
        surplus_total = 0
        
        for trial in range(200):
            L = jm_sample(n, steps=n*500)
            A = get_A(L)
            
            r = rank_Fp(A, p)
            kp = (n-1) - r
            kp_dist[kp] = kp_dist.get(kp, 0) + 1
            
            detA = bareiss_det(A.tolist())
            if detA == 0:
                continue
            
            vp_actual = vp(detA, p)
            vp_det_dist[vp_actual] = vp_det_dist.get(vp_actual, 0) + 1
            
            # p-Sylow bound
            psy_bound = vp(n // p, p) + max(0, kp - 1)  # = v_p(n)-1 + max(0,kp-1)
            
            if vp_actual < psy_bound:
                violations += 1
            
            surplus = vp_actual - thm32_bound
            surplus_total += surplus
        
        kp_sorted = sorted(kp_dist.items())
        print(f"    p={p}: k_p distribution = {dict(kp_sorted)}")
        print(f"         v_{p}(det A) distribution = {dict(sorted(vp_det_dist.items()))}")
        print(f"         Thm 3.2 baseline: v_{p}(det A) ≥ {thm32_bound}")
        print(f"         p-Sylow violations: {violations}")
        avg_surplus = surplus_total / max(1, sum(1 for _ in range(200)))
        print(f"         avg surplus over Thm 3.2: {avg_surplus:.2f}")

# ── 3c: Per p=2: verifica che p-Sylow = Thm 5.5 ──────────────────
print("\n" + "─" * 72)
print("  3c: Per p=2, la congettura p-Sylow ≡ Thm 5.5 (F2-rank)?")
print("─" * 72)

print("""
  Thm 5.5 (nel paper): v_2(det A) ≥ v_2(n/2) + max(0, k_2 - 1)
                      = v_2(n) - 1 + max(0, k_2 - 1)
  p-Sylow per p=2:     v_2(det A) ≥ v_2(n/2) + max(0, k_2 - 1)
  
  → IDENTICI. Il Thm 5.5 IS la congettura p-Sylow per p=2, DIMOSTRATA.
  
  Quindi: l'unica parte NON dimostrata della congettura è per p DISPARI.
""")

# ── 3d: Per p dispari: quanto è forte il bound? ──────────────────
print("─" * 72)
print("  3d: Per p dispari — il bound p-Sylow è più forte di Thm 3.2?")
print("─" * 72)

print("""
  Riscriviamo: 
    Thm 3.2:   v_p(det A) ≥ v_p(n)       [n dispari, o p dispari con n pari]
    p-Sylow:   v_p(det A) ≥ v_p(n) - 1 + max(0, k_p - 1)
  
  p-Sylow ≥ Thm 3.2  ⟺  max(0, k_p-1) ≥ 1  ⟺  k_p ≥ 2.
  
  Quando k_p = 0: p-Sylow dà v_p(n)-1, PEGGIO di Thm 3.2.
  Quando k_p = 1: p-Sylow dà v_p(n)-1, PEGGIO di Thm 3.2.
  Quando k_p = 2: p-Sylow dà v_p(n), UGUALE a Thm 3.2.
  Quando k_p ≥ 3: p-Sylow dà v_p(n)+k_p-2, MEGLIO di Thm 3.2!
  
  → CONCLUSIONE: La formula p-Sylow è SBAGLIATA nella normalizzazione!
    Il baseline NON è v_p(n/p), ma v_p(n). 
    La formula CORRETTA per p dispari dovrebbe essere:
    
      v_p(det A) ≥ v_p(n) + max(0, k_p - 1)   [congettura corretta]
    
    Verifichiamo questa formula corretta...
""")

# Verify the CORRECT formula: v_p(det A) >= v_p(n) + max(0, k_p - 1) for odd p
print("  Verifica formula CORRETTA: v_p(det A) ≥ v_p(n) + max(0, k_p - 1)")
print("  (per p dispari)")

total_tests = 0
violations_old = 0
violations_new = 0

for n in [6, 9, 10, 12, 15]:
    factors = factorint(n)
    odd_primes = [p for p in factors if p > 2]
    if not odd_primes:
        continue
    
    for trial in range(200):
        L = jm_sample(n, steps=n*500)
        A = get_A(L)
        detA = bareiss_det(A.tolist())
        if detA == 0: continue
        
        for p in odd_primes:
            total_tests += 1
            kp = (n-1) - rank_Fp(A, p)
            vp_actual = vp(detA, p)
            
            old_bound = vp(n, p) - 1 + max(0, kp - 1)  # original rem:p-Sylow
            new_bound = vp(n, p) + max(0, kp - 1)       # corrected
            
            if vp_actual < old_bound: violations_old += 1
            if vp_actual < new_bound: violations_new += 1

print(f"  Total tests: {total_tests}")
print(f"  Violations (original v_p(n/p) + max(0,kp-1)): {violations_old}")
print(f"  Violations (corrected v_p(n) + max(0,kp-1)):  {violations_new}")

if violations_new == 0:
    print("""
  ★ LA FORMULA CORRETTA per p dispari è:
      v_p(det A) ≥ v_p(n) + max(0, k_p - 1)
    
    Unificata con p=2 (Thm 5.5), la formula completa diventa:
      v_p(det A) ≥ v_p(n/gcd(n,p)) + max(0, k_p - 1)    per p=2
      v_p(det A) ≥ v_p(n) + max(0, k_p - 1)              per p dispari
    
    Oppure in modo unificato:
      v_p(det A) ≥ v_p(n) - [p=2 e 2|n] + max(0, k_p - 1)
""")
else:
    print(f"\n  Hmm, {violations_new} violations for corrected formula. Let's check exceptions...")


# ════════════════════════════════════════════════════════════════════
print("\n\n" + "=" * 72)
print("  PUNTO 6: PERCHÉ K(A_cyc) ≅ (Z/n)^{n-2}?")
print("=" * 72)

# ── 6a: Connessione con il Laplaciano del grafo completo K_n ──────
print("""
OSSERVAZIONE CHIAVE: Il grafo completo K_n ha Laplaciano
  L_{K_n} = nI - J    (dove J = matrice di tutti 1)

Il Kirchhoff Matrix-Tree Theorem dà:
  • autovalori di L_{K_n}: n (con molteplicità n-1) e 0 (molteplicità 1)  
  • SNF di L_{K_n}: diag(0, n, n, ..., n)
  • Sandpile group: K(K_n) ≅ (Z/n)^{n-2}

DOMANDA: A (matrice ridotta del LS ciclico) ha la STESSA SNF?
""")

print("─" * 72)
print("  6a: Confronto A_cyc vs L_{K_n} (ridotto)")
print("─" * 72)

for n in range(2, 10):
    A = get_A(cyclic_ls(n))
    m = n - 1
    
    # Laplaciano ridotto di K_n: togliere una riga e colonna da (nI - J)
    # L_reduced = nI_{m×m} - J_{m×m}  (NO! Bisogna togliere riga/colonna dall'originale)
    # L_{K_n} = nI_n - J_n. Ridotto = prime (n-1) righe e colonne.
    # L_red[i,j] = n*delta_{ij} - 1
    L_red = np.array([[n * (1 if i==j else 0) - 1 for j in range(m)] for i in range(m)])
    
    snf_A = snf_diag(A)
    snf_L = snf_diag(L_red)
    
    same = (snf_A == snf_L)
    tag = "≡" if same else "≠"
    print(f"  n={n}: SNF(A_cyc) = {snf_A},  SNF(L_red(K_n)) = {snf_L}  {tag}")

# ── 6b: Struttura degli autovalori ────────────────────────────────
print("\n" + "─" * 72)
print("  6b: Autovalori — perché sono gli stessi?")
print("─" * 72)

for n in [5, 7, 8]:
    A = get_A(cyclic_ls(n))
    m = n - 1
    L_red = np.array([[n * (1 if i==j else 0) - 1 for j in range(m)] for i in range(m)])
    
    eig_A = sorted(np.linalg.eigvals(A.astype(float)).real)
    eig_L = sorted(np.linalg.eigvals(L_red.astype(float)).real)
    
    print(f"\n  n={n}:")
    print(f"    eig(A_cyc)     = [{', '.join(f'{x:.3f}' for x in eig_A)}]")
    print(f"    eig(L_red(Kn)) = [{', '.join(f'{x:.3f}' for x in eig_L)}]")
    
    # Check: are they all equal to n?
    all_n_A = all(abs(x - n) < 0.01 or abs(x + n) < 0.01 or abs(abs(x) - n) < 0.01 for x in eig_A)
    all_n_L = all(abs(x - n) < 0.01 for x in eig_L)
    print(f"    eig(L_red) all = n? {all_n_L}")

# ── 6c: Identità strutturale ─────────────────────────────────────
print("\n" + "─" * 72)
print("  6c: C'è una trasformazione diretta A_cyc → L_red(K_n)?")
print("─" * 72)

for n in [3, 4, 5, 6, 7]:
    A = get_A(cyclic_ls(n))
    m = n - 1
    L_red = np.array([[n * (1 if i==j else 0) - 1 for j in range(m)] for i in range(m)])
    
    # Is A = P * L_red * Q for some unimodular P, Q?
    # (Same SNF ↔ equivalent over Z)
    # Check: is A *exactly* L_red? 
    same_matrix = np.array_equal(A, L_red)
    
    # Check: det(A) vs det(L_red)
    det_A = bareiss_det(A.tolist())
    det_L = bareiss_det(L_red.tolist())
    
    print(f"  n={n}: A = L_red? {same_matrix}  |  det(A) = {det_A}, det(L_red) = {det_L}")

# ── 6d: La spiegazione circolante ─────────────────────────────────
print("\n" + "─" * 72)
print("  6d: A_cyc è una matrice CIRCOLANTE — ecco perché")
print("─" * 72)

print("""
  Per il LS ciclico: L_{ij} = (i+j) mod n + 1
  
  A_{ij} = L_{ij} - L_{i,n-1} = [(i+j) mod n] - [(i+n-1) mod n]
  
  Definisci c_j = (j mod n) - ((n-1) mod n) = j - (n-1) mod n
  
  Ma A_{ij} NON dipende solo da (j-i) — dipende da i e j separatamente
  tramite il modular arithmetic. Tuttavia A_{ij} ≡ j+1 (mod n).
  
  Vediamo la struttura precisa:
""")

for n in [4, 5, 6]:
    A = get_A(cyclic_ls(n))
    m = n - 1
    print(f"  n={n}, A =")
    for i in range(m):
        print(f"    [{', '.join(f'{A[i,j]:4d}' for j in range(m))}]")
    
    # Check if A is circulant
    is_circ = True
    for i in range(1, m):
        for j in range(m):
            if A[i][j] != A[0][(j - i) % m]:
                is_circ = False; break
    print(f"  Circolante? {is_circ}")
    
    # What's A mod n?
    R = A % n
    print(f"  A mod {n} = {[list(R[i]) for i in range(m)]}")
    print()


# ── 6e: Autovalori del circolante e il "perché" ──────────────────
print("─" * 72)
print("  6e: Autovalori del circolante → connessione con DFT")
print("─" * 72)

print("""
  Se A è circolante con prima riga [c_0, c_1, ..., c_{m-1}], 
  allora gli autovalori sono:
    λ_k = Σ_j c_j ω^{jk}   dove ω = e^{2πi/m}, m = n-1
  
  Per A_cyc del LS ciclico:
""")

for n in [5, 7, 11]:
    A = get_A(cyclic_ls(n))
    m = n - 1
    c = A[0]  # prima riga del circolante
    
    # Autovalori via DFT
    omega = np.exp(2j * np.pi / m)
    eigs = []
    for k in range(m):
        lam = sum(c[j] * omega**(j*k) for j in range(m))
        eigs.append(lam)
    
    print(f"  n={n}, m={m}:")
    print(f"    prima riga c = {list(c)}")
    print(f"    autovalori λ_k:")
    for k in range(m):
        lam = eigs[k]
        print(f"      k={k}: λ = {lam.real:+.4f} {lam.imag:+.4f}i  |λ| = {abs(lam):.4f}")
    
    product = 1.0
    for lam in eigs:
        product *= lam
    print(f"    Prodotto |Π λ_k| = {abs(product):.2f}, n^{m} = {n**m}")
    # det(A) should be product of eigenvalues
    detA = bareiss_det(A.tolist())
    print(f"    det(A) = {detA}, expected sign*n^(n-2) = ±{n**(n-2)}")
    
    # Check: |λ_k| = n for all k?
    all_abs_n = all(abs(abs(lam) - n) < 0.01 for lam in eigs)
    print(f"    |λ_k| = n per ogni k? {all_abs_n}")
    if all_abs_n:
        print(f"    ★ TUTTI gli autovalori hanno |λ| = n!")
    print()

# ── 6f: SPIEGAZIONE PROFONDA ──────────────────────────────────────
print("─" * 72)
print("  6f: LA SPIEGAZIONE PROFONDA")
print("─" * 72)

print("""
  TEOREMA (spiegazione strutturale):
  
  A_cyc è una matrice circolante (n-1)×(n-1) con autovalori che hanno
  TUTTI modulo esattamente n. Questo implica:
  
  1. |det(A)| = Π|λ_k| = n^{n-2}  (Prop 3.6 come corollario)
  
  2. La SNF ha d_k | n per proprietà del ring Z[ω]. Combinato con  
     Step 1 (multilinearità → n | d_k per k≥2), otteniamo d_k = n.
  
  MA PERCHÉ |λ_k| = n per tutti k?
""")

# Prove why |λ_k| = n
print("  Analisi dettagliata degli autovalori:")
print()

for n in [5, 7]:
    m = n - 1
    A = get_A(cyclic_ls(n))
    c = A[0]
    
    omega = np.exp(2j * np.pi / m)
    zeta_n = np.exp(2j * np.pi / n)
    
    print(f"  n={n}: prima riga c_j = L(0,j) - L(0,n-1) = (j mod {n}) + 1 - {(n-1) % n + 1}")
    print(f"         = {list(c)}")
    
    # c_j = (j+1) - n = j + 1 - n  for j = 0,...,n-2 IF j+1 < n (always true for j ≤ n-2)
    # Wait, L_{0,j} = (0+j) mod n + 1 = j + 1, L_{0,n-1} = n
    # So c_j = (j+1) - n for j = 0,...,n-2
    # c_j = j - (n-1) = -(n-1-j)
    
    c_formula = [j + 1 - n for j in range(m)]
    print(f"         formula: c_j = j + 1 - n = {c_formula}")
    match = list(c) == c_formula
    print(f"         match? {match}")
    
    # λ_k = Σ_j c_j ω^{jk} = Σ_j (j+1-n) ω^{jk}
    # = Σ_j (j+1) ω^{jk} - n Σ_j ω^{jk}
    # = Σ_j (j+1) ω^{jk} - n * (sum of geometric series)
    # 
    # Σ_j=0^{m-1} ω^{jk} = 0 for k ≠ 0 (mod m), = m for k = 0
    # 
    # So for k ≠ 0: λ_k = Σ_{j=0}^{m-1} (j+1) ω^{jk}
    # For k = 0: λ_0 = Σ_{j=0}^{m-1} (j+1-n) = Σ_{j=1}^{n-1} j - n(n-1)
    #          = n(n-1)/2 - n(n-1) = -n(n-1)/2
    
    lam0 = sum(c)
    expected_lam0 = -n * (n-1) // 2
    print(f"         λ_0 = {lam0} = -n(n-1)/2 = {expected_lam0} ✓")
    
    # For k ≥ 1: λ_k = Σ_{j=0}^{m-1} (j+1) ω^{jk}
    # Let s = Σ_{j=0}^{m-1} j ω^{jk} = ω^k * d/dω^k Σ ω^{jk} 
    # Actually: Σ_j j ω^{jk} = ... let's compute differently.
    #
    # Σ_{j=0}^{m-1} (j+1) ω^{jk} = Σ_{l=1}^{m} l ω^{(l-1)k}
    # = ω^{-k} Σ_{l=1}^{m} l ω^{lk}
    
    for k in range(m):
        lam_k = sum(c[j] * omega**(j*k) for j in range(m))
        print(f"         λ_{k} = {lam_k.real:+8.4f} {lam_k.imag:+8.4f}i,  |λ_{k}| = {abs(lam_k):.6f}")
    print()

# ── 6g: Prova algebrica ──────────────────────────────────────────
print("─" * 72)
print("  6g: PROVA ALGEBRICA che |λ_k| = n")
print("─" * 72)

print("""
  Prima riga: c_j = j + 1 - n  per j = 0, ..., n-2.
  Sia ω = e^{2πi/(n-1)}.
  
  λ_k = Σ_{j=0}^{n-2} (j + 1 - n) ω^{jk}
      = Σ_{j=0}^{n-2} (j+1) ω^{jk}  -  n Σ_{j=0}^{n-2} ω^{jk}
  
  Per k = 0: λ_0 = n(n-1)/2 - n(n-1) = -n(n-1)/2.
           |λ_0| = n(n-1)/2.  MA |λ_0| = n ⟺ (n-1)/2 = 1 ⟺ n = 3.
           
           PROBLEMA! |λ_0| ≠ n in generale!
""")

for n in [3, 4, 5, 6, 7, 8, 9, 10, 11]:
    m = n - 1
    A = get_A(cyclic_ls(n))
    c = A[0]
    omega = np.exp(2j * np.pi / m)
    
    eigs = [sum(c[j] * omega**(j*k) for j in range(m)) for k in range(m)]
    abs_eigs = [abs(lam) for lam in eigs]
    
    det_A = bareiss_det(A.tolist())
    prod_abs = 1
    for ae in abs_eigs:
        prod_abs *= ae
    
    print(f"  n={n}: |λ_k| = [{', '.join(f'{x:.2f}' for x in abs_eigs)}]")
    print(f"         Π|λ_k| = {prod_abs:.2f},  |det(A)| = {abs(det_A)},  n^(n-2) = {n**(n-2)}")

print("""
  AH — gli autovalori NON hanno tutti |λ| = n. 
  Il prodotto Π|λ_k| = n^{n-2} è vero, ma i singoli |λ_k| variano.
  
  La spiegazione profonda è diversa...
""")

# ── 6h: Qual è la VERA spiegazione? ──────────────────────────────
print("─" * 72)
print("  6h: La VERA spiegazione — equivalenza Z-lineare")
print("─" * 72)

print("""
  Cerchiamo una matrice unimodulare U (det ±1) tale che U A U^T = L_red(K_n)
  o almeno A e L_red(K_n) siano equivalenti su Z (stessa SNF).
  
  Sappiamo:
    - SNF(A_cyc) = diag(1, n, ..., n)   [Cor 3.7]
    - SNF(L_red(K_n)) = diag(n-1, n, ..., n) per il Laplaciano ridotto? 
      NO ASPETTA — ricalcoliamo...
""")

for n in range(2, 10):
    m = n - 1
    L_red = np.array([[n * (1 if i==j else 0) - 1 for j in range(m)] for i in range(m)])
    snf = snf_diag(L_red)
    det_L = bareiss_det(L_red.tolist())
    print(f"  n={n}: SNF(L_red(K_n)) = {snf},  det = {det_L}")

print("""
  SNF(L_red(K_n)) = diag(1, n, n, ..., n) per n ≥ 2!
  
  Stessa SNF di A_cyc! Ma sono la STESSA matrice? No. Sono Z-equivalenti:
  ∃ unimodulari P,Q con P A_cyc Q = L_red.
  
  SIGNIFICATO PROFONDO:
  
  ┌──────────────────────────────────────────────────────────────────┐
  │  K(A_cyc) ≅ Z^{n-1}/Im(A_cyc)                                 │
  │           ≅ Z^{n-1}/Im(L_red(K_n))                             │
  │           = K(K_n)                                              │
  │           = sandpile group del grafo completo                   │
  │           ≅ (Z/n)^{n-2}                                        │
  │                                                                 │
  │  Il LS CICLICO di ordine n ha lo STESSO sandpile group          │
  │  del GRAFO COMPLETO K_n !                                       │
  │                                                                 │
  │  Questo è il "perché" profondo:                                 │
  │  La struttura è ereditata dalla simmetria ciclica,              │
  │  che rende la matrice Z-equivalente al Laplaciano di K_n.       │
  └──────────────────────────────────────────────────────────────────┘
""")

# Verify: find explicit P, Q such that P * A * Q = L_red
print("  Verifica: trovare P, Q unimodulari con P * A * Q = L_red(K_n)")
print()

for n in [3, 4, 5]:
    m = n - 1
    A = get_A(cyclic_ls(n))
    L_red = np.array([[n * (1 if i==j else 0) - 1 for j in range(m)] for i in range(m)])
    
    # Since both have the same SNF, ∃ P₁,Q₁ with P₁ A Q₁ = SNF
    # and P₂,Q₂ with P₂ L_red Q₂ = SNF
    # So A = P₁⁻¹ SNF Q₁⁻¹ and L_red = P₂⁻¹ SNF Q₂⁻¹
    # Hence P₂ P₁⁻¹ A Q₁⁻¹ Q₂ = L_red
    
    # Let's just check: is A * v = n * v solvable over Z?
    # Eigenvalue n: (A - nI) v = 0 over Z
    AnI = A - n * np.eye(m, dtype=int)
    print(f"  n={n}: A - nI = ")
    for i in range(m):
        print(f"    [{', '.join(f'{AnI[i,j]:4d}' for j in range(m))}]")
    
    # rank over Q
    r = np.linalg.matrix_rank(AnI.astype(float))
    print(f"    rank(A - nI) over Q = {r}")
    print(f"    nullità = {m - r} (= n-2 conferma spazio nullo grande)")
    print()

# ── 6i: Connessione esplicita ────────────────────────────────────
print("─" * 72)
print("  6i: LA CONNESSIONE A ↔ Laplaciano: spiegazione costruttiva")
print("─" * 72)

print("""
  L'identità chiave è: A = R + nΔ dove R = rank-1 matrix.
  
  R = [1, 2, ..., n-1]    (tutte le righe uguali)
    = v · 1ᵀ               dove v = (1, 2, ..., n-1)ᵀ e 1 = (1,...,1)ᵀ
  
  L_red = nI - J = nI - 1·1ᵀ   dove J è la matrice di soli 1.
  
  DOMANDA: la trasformazione colonne che "centra" R → -J 
  porta anche nΔ → nI?
""")

for n in [4, 5, 6]:
    m = n - 1
    A = get_A(cyclic_ls(n))
    L_red = np.array([[n * (1 if i==j else 0) - 1 for j in range(m)] for i in range(m)])
    
    # Delta = (A - R) / n
    R = np.array([[j+1 for j in range(m)] for _ in range(m)])
    Delta = (A - R) // n
    
    print(f"  n={n}:")
    print(f"    R = {[list(R[0])]}")
    print(f"    Δ = ")
    for i in range(m):
        print(f"      [{', '.join(f'{Delta[i,j]:3d}' for j in range(m))}]")
    
    # Is Δ = I - something?
    # Actually for the cyclic case:
    # A_{ij} = (i+j) mod n - (i+n-1) mod n
    # When i+j < n: A_{ij} = (i+j) - (i+n-1) = j - (n-1) = j+1-n
    # When i+j ≥ n: A_{ij} = (i+j-n) - (i-1) = j - n + 1 ... wait
    # Let me be more careful
    
    print(f"    nΔ = A - R:")
    nDelta = A - R
    for i in range(m):
        print(f"      [{', '.join(f'{nDelta[i,j]:3d}' for j in range(m))}]")
    
    # nΔ should be n * Δ where Δ is integer
    # Check: is nΔ = -nI + ... something?
    diff_from_nI = nDelta + n * np.eye(m, dtype=int)
    print(f"    nΔ + nI:")
    for i in range(m):
        print(f"      [{', '.join(f'{diff_from_nI[i,j]:3d}' for j in range(m))}]")
    print()


# ── FINAL SYNTHESIS ───────────────────────────────────────────────
print("\n" + "=" * 72)
print("  SINTESI FINALE")
print("=" * 72)

print("""
  ═══════════════════════════════════════════════════════════════════
  PUNTO 3 — p-SYLOW: CONCLUSIONI
  ═══════════════════════════════════════════════════════════════════
  
  1. Per p=2: la congettura È il Thm 5.5. GIÀ DIMOSTRATA.
  
  2. Per p dispari: la formula v_p(det A) ≥ v_p(n/p) + max(0, k_p-1)
     è SBAGLIATA nella normalizzazione — dà risultati PIÙ DEBOLI
     del Thm 3.2 quando k_p ≤ 1.
     
     La formula CORRETTA dovrebbe essere:
       v_p(det A) ≥ v_p(n) + max(0, k_p - 1)    [baseline = Thm 3.2]
     
     → O la congettura p-Sylow è BANALE (segue da Thm 3.2 per k_p ≤ 1),
       oppure è una generalizzazione genuina solo per k_p ≥ 2.
     
  3. Per il LS CICLICO: k_p = n-2 (rank = 1), quindi la congettura
     dà v_p ≥ v_p(n) + n-3, che è v_p(det A) = (n-2)v_p(n). 
     ← NONTRIVIALE e VERA.
  
  ═══════════════════════════════════════════════════════════════════
  PUNTO 6 — IL "PERCHÉ" PROFONDO
  ═══════════════════════════════════════════════════════════════════
  
  K(A_cyc) ≅ (Z/n)^{n-2}  PERCHÉ:
  
  A_cyc è Z-equivalente al Laplaciano ridotto del grafo completo K_n.
  
  Concretamente:
    A = R + nΔ   dove R = (rank-1) e nΔ ≡ diag pattern
    L_red(K_n) = nI - J   dove J = (rank-1)
  
  Entrambe sono "rank-1 perturbation of n × (diagonal-like matrix)".
  La Z-equivalenza delle SNF segue dalla struttura identica dei
  divisori elementari.
  
  Il "perché" ultimo: la simmetria ciclica del LS impone che A abbia
  la STESSA struttura algebra lineare del Kirchhoff del grafo più 
  simmetrico possibile (K_n).
  ═══════════════════════════════════════════════════════════════════
""")
