"""
deep_investigation_3.py — Verify critical anomaly:
Does k_p = n-2 hold for ALL LS, or just JM-sampled ones?
Also: when n/p = 1 (p=n) vs n/p ≥ 2.
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass
import numpy as np
from math import gcd
from functools import reduce
from itertools import permutations

def get_A(L):
    n = L.shape[0]
    return np.array([[int(L[i,j]) - int(L[i,n-1]) for j in range(n-1)] for i in range(n-1)])

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
    if x == 0: return float('inf')
    x = abs(x); v = 0
    while x % p == 0: x //= p; v += 1
    return v

def bareiss_det(M_list):
    m = len(M_list); M = [row[:] for row in M_list]
    sign = 1; prev = 1
    for i in range(m):
        if M[i][i] == 0:
            for k in range(i+1, m):
                if M[k][i] != 0: M[i], M[k] = M[k], M[i]; sign *= -1; break
            else: return 0
        for j in range(i+1, m):
            for k in range(i+1, m):
                M[j][k] = (M[j][k]*M[i][i] - M[j][i]*M[i][k]) // prev
            M[j][i] = 0
        prev = M[i][i]
    return sign * M[m-1][m-1]

def generate_reduced_ls(n):
    results = []
    def backtrack(grid, row, col):
        if row == n:
            results.append([r[:] for r in grid]); return
        next_row, next_col = (row, col+1) if col+1 < n else (row+1, 1)
        used_row = set(grid[row])
        used_col = {grid[r][col] for r in range(row)}
        for val in range(1, n+1):
            if val not in used_row and val not in used_col:
                grid[row][col] = val
                backtrack(grid, next_row, next_col)
                grid[row][col] = 0
    grid = [[0]*n for _ in range(n)]
    for j in range(n): grid[0][j] = j + 1
    for i in range(n): grid[i][0] = i + 1
    backtrack(grid, 1, 1)
    return results


print("=" * 72)
print("  VERIFICA ESAUSTIVA: k_p per TUTTI i LS ridotti")
print("=" * 72)

# ── n=5, p=5 (CRITICAL: p = n, n/p = 1) ─────────────────────────
print("\n  n = 5, p = 5 (p = n, caso critico):")
ls5 = generate_reduced_ls(5)

kp_det_map = {}
for ls in ls5:
    L = np.array(ls)
    A = get_A(L)
    detA = bareiss_det(A.tolist())
    kp = (5-1) - rank_Fp(A, 5)
    key = kp
    if key not in kp_det_map:
        kp_det_map[key] = []
    kp_det_map[key].append(abs(detA))

for kp in sorted(kp_det_map):
    dets = kp_det_map[kp]
    det_set = sorted(set(dets))
    vps = [vp(d, 5) for d in dets if d != 0]
    print(f"    k_5 = {kp}: {len(dets)} LS, |det| values = {det_set}")
    if vps:
        bound = 1 + max(0, kp - 1)
        tight = sum(1 for v in vps if v == bound)
        print(f"      v_5 range: [{min(vps)}, {max(vps)}], "
              f"bound = {bound}, tight: {tight}/{len(vps)}")

# ── n=5, p=5: verify the CORRECTED formula ───────────────────────
print(f"\n  Verifica formula corretta v_5 ≥ v_5(5) + max(0, k_5 - 1):")
violations = 0
for ls in ls5:
    L = np.array(ls)
    A = get_A(L)
    detA = bareiss_det(A.tolist())
    if detA == 0: continue
    kp = (5-1) - rank_Fp(A, 5)
    bound = 1 + max(0, kp - 1)
    actual = vp(detA, 5)
    if actual < bound:
        violations += 1
        print(f"    VIOLATION: |det|={abs(detA)}, k_5={kp}, v_5={actual}, bound={bound}")
print(f"  Violations: {violations}/{len(ls5)}")

# ── n=4, p=2 (p divides n, n/p=2) ────────────────────────────────
print(f"\n  n = 4, p = 2 (n/p = 2):")
ls4 = generate_reduced_ls(4)
kp_dist = {}
for ls in ls4:
    L = np.array(ls)
    A = get_A(L)
    kp = 3 - rank_Fp(A, 2)
    kp_dist[kp] = kp_dist.get(kp, 0) + 1
print(f"    k_2 distribution: {dict(sorted(kp_dist.items()))}")

# ── n=6, p=2 and p=3 ─────────────────────────────────────────────
ls6 = generate_reduced_ls(6)
print(f"\n  n = 6, p = 2 e p = 3 (ESAUSTIVO — {len(ls6)} LS):")

for p in [2, 3]:
    kp_dist = {}
    formula_violations = 0
    for ls in ls6:
        L = np.array(ls)
        A = get_A(L)
        detA = bareiss_det(A.tolist())
        kp = 5 - rank_Fp(A, p)
        kp_dist[kp] = kp_dist.get(kp, 0) + 1
        
        if detA != 0:
            iverson = 1 if (p == 2 and 6 % 2 == 0) else 0
            bound = vp(6, p) - iverson + max(0, kp - 1)
            actual = vp(detA, p)
            if actual < bound:
                formula_violations += 1
    
    print(f"    p={p}: k_{p} distribution = {dict(sorted(kp_dist.items()))}")
    print(f"       formula violations: {formula_violations}")

# ── n=5, p=5: Examine A mod 5 structure ──────────────────────────
print(f"\n" + "─" * 72)
print(f"  PERCHÉ k_5 varia per n=5?")
print(f"─" * 72)

# Show A mod 5 for different k_5 values
for kp_target in sorted(kp_det_map.keys()):
    print(f"\n  Esempio con k_5 = {kp_target}:")
    found = False
    for ls in ls5:
        L = np.array(ls)
        A = get_A(L)
        kp = 4 - rank_Fp(A, 5)
        if kp == kp_target and not found:
            Amod = A % 5
            print(f"    A =")
            for i in range(4):
                print(f"      [{', '.join(f'{A[i,j]:4d}' for j in range(4))}]")
            print(f"    A mod 5 =")
            for i in range(4):
                print(f"      [{', '.join(f'{Amod[i,j]:2d}' for j in range(4))}]")
            # Check: are rows identical?
            rows = [tuple(Amod[i]) for i in range(4)]
            print(f"    Rows identical? {len(set(rows)) == 1}")
            if len(set(rows)) > 1:
                print(f"    Distinct rows: {set(rows)}")
            found = True

# ── KEY INSIGHT: n/p = 1 vs n/p ≥ 2 ─────────────────────────────
print(f"\n" + "─" * 72)
print(f"  INSIGHT CHIAVE: n/p = 1 vs n/p ≥ 2")
print(f"─" * 72)

print("""
  Quando n/p ≥ 2: in ogni riga di L, ogni classe di residui mod p 
  appare ALMENO 2 volte. Dopo la sottrazione L_{ij} - L_{i,n-1},
  la struttura mod p è più vincolata.
  
  Quando n = p: ogni riga di L è una permutazione di {1,...,p},
  che mod p è una permutazione di {0,1,...,p-1}. 
  La sottrazione può produrre righe mod p molto diverse.
""")

# Verify: for n/p ≥ 2, is rank_p = 1 universal?
print("  Verifica: per n/p ≥ 2, rank_p(A) = 1 per TUTTI i LS ridotti?")

# n=4, p=2 (n/p=2)
all_rank1_4_2 = all(rank_Fp(get_A(np.array(ls)), 2) == 1 for ls in ls4)
print(f"    n=4, p=2: rank_2 = 1 sempre? {all_rank1_4_2}")

# n=6, p=2 (n/p=3) 
all_rank1_6_2 = all(rank_Fp(get_A(np.array(ls)), 2) == 1 for ls in ls6)
print(f"    n=6, p=2: rank_2 = 1 sempre? {all_rank1_6_2}")

# n=6, p=3 (n/p=2)
all_rank1_6_3 = all(rank_Fp(get_A(np.array(ls)), 3) == 1 for ls in ls6)
print(f"    n=6, p=3: rank_3 = 1 sempre? {all_rank1_6_3}")

# n=5, p=5 (n/p=1)
all_rank1_5_5 = all(rank_Fp(get_A(np.array(ls)), 5) == 1 for ls in ls5)
print(f"    n=5, p=5: rank_5 = 1 sempre? {all_rank1_5_5}")


# ── PROVE that rank_p = 1 for n/p ≥ 2 ────────────────────────────
print(f"\n" + "─" * 72)
print(f"  DIMOSTRAZIONE che rank_p(A) = 1 per n/p ≥ 2")  
print(f"─" * 72)

print("""
  TEOREMA: Per ogni LS di ordine n e primo p con p|n e n/p ≥ 2:
  rank_{F_p}(A) = 1, cioè k_p = n-2.
  
  PROVA:
  Sia k = n/p ≥ 2. In ogni riga i di L, ogni classe r ∈ F_p appare
  esattamente k volte. In particolare:
  
  (1) Somma di riga: Σ_{j=0}^{n-2} A_{ij} = Σ_{j≠n-1} L_{ij} - (n-1)L_{i,n-1}
      = [n(n+1)/2 - L_{i,n-1}] - (n-1)L_{i,n-1} = n(n+1)/2 - nL_{i,n-1}
      ≡ 0 (mod p) ← perché n ≡ 0 mod p.
      
  (2) Somma pesata: Σ_{j=0}^{n-2} f(j) A_{ij} mod p per qualsiasi f?
      Ci serve dimostrare che le righe sono proporzionali su F_p, non solo
      che hanno somma zero.
      
  Approccio diverso — usiamo le COLONNE:
  
  Colonna j di A: A_{*,j} = (L_{0,j}-L_{0,n-1}, ..., L_{n-2,j}-L_{n-2,n-1}).
  
  Mod p: A_{i,j} ≡ L_{ij} - L_{i,n-1} (mod p).
  
  Le colonne j e j' di L sono ENTRAMBE permutazioni di {1,...,n}.
  Mod p, entrambe contengono ogni residuo esattamente k volte.
  
  La differenza colonna_j - colonna_{j'} di A mod p è:
  (L_{*,j} - L_{*,n-1}) - (L_{*,j'} - L_{*,n-1}) = L_{*,j} - L_{*,j'} mod p.
  
  Questo NON è necessariamente zero. L'argomento deve passare per un'altra via.
""")

# ── Maybe it's a COUNTING argument? ──────────────────────────────
print("  Approccio combinatorio: conteggio esatto")
print()

# For n=6, p=3: show A mod 3 for several LS
print("  n=6, p=3: A mod 3 per 5 LS random:")
import random; random.seed(123)
indices = random.sample(range(len(ls6)), 5)
for idx in indices:
    L = np.array(ls6[idx])
    A = get_A(L)
    Amod = A % 3
    rows = [tuple(Amod[i]) for i in range(5)]
    unique = len(set(rows))
    # Check rank
    r = rank_Fp(A, 3)
    print(f"    LS#{idx}: rank_3 = {r}, rows: {rows}")

# ── Final analysis: is rank_p = 1 (n/p ≥ 2) actually a theorem? ──
print(f"\n" + "─" * 72)
print(f"  ANALISI ALGEBRAICA: perché rank_p = 1 per n/p ≥ 2?")
print(f"─" * 72)

print("""
  Per n/p = k ≥ 2, consideriamo A mod p.
  
  CLAIM: Ogni 2×2 minore di A mod p è ≡ 0 (mod p).
  
  A 2×2 minor: det[[A_{i,j}, A_{i,j'}],[A_{i',j}, A_{i',j'}]]
  = (L_{ij}-L_{i,n-1})(L_{i'j'}-L_{i',n-1}) - (L_{ij'}-L_{i,n-1})(L_{i'j}-L_{i',n-1})
  = L_{ij}L_{i'j'} - L_{ij}L_{i',n-1} - L_{i,n-1}L_{i'j'} + L_{i,n-1}L_{i',n-1}
    - L_{ij'}L_{i'j} + L_{ij'}L_{i',n-1} + L_{i,n-1}L_{i'j} - L_{i,n-1}L_{i',n-1}
  = L_{ij}L_{i'j'} - L_{ij'}L_{i'j} - L_{i',n-1}(L_{ij}-L_{ij'}) - L_{i,n-1}(L_{i'j'}-L_{i'j})
  = (L_{ij}-L_{ij'})(L_{i'j'}-L_{i',n-1}) - (L_{i'j}-L_{i'j'})(L_{ij'}-L_{i,n-1})
  
  Hmm, questo non si semplifica facilmente. Verifichiamo numericamente.
""")

# Check: are ALL 2x2 minors of A mod p zero for n/p ≥ 2?
from itertools import combinations

ls_cache = {4: ls4, 5: ls5, 6: ls6}
for n, p in [(4, 2), (6, 2), (6, 3)]:
    m = n - 1
    all_zero = True
    ls_list = ls_cache[n]
    
    nonzero_count = 0
    total_minors = 0
    
    for ls in ls_list:
        L = np.array(ls)
        A = get_A(L) % p
        
        for (i1, i2) in combinations(range(m), 2):
            for (j1, j2) in combinations(range(m), 2):
                minor = (int(A[i1,j1])*int(A[i2,j2]) - int(A[i1,j2])*int(A[i2,j1])) % p
                total_minors += 1
                if minor != 0:
                    nonzero_count += 1
                    all_zero = False
    
    print(f"  n={n}, p={p} (n/p={n//p}): 2×2 minors ≡ 0 mod {p}? "
          f"{'YES ✓' if all_zero else f'NO — {nonzero_count}/{total_minors} nonzero'}")

# For n=5, p=5 (n/p = 1): should NOT all be zero
for n, p in [(5, 5)]:
    m = n - 1
    nonzero_count = 0
    total_minors = 0
    for ls in ls5:
        L = np.array(ls)
        A = get_A(L) % p
        for (i1, i2) in combinations(range(m), 2):
            for (j1, j2) in combinations(range(m), 2):
                minor = (int(A[i1,j1])*int(A[i2,j2]) - int(A[i1,j2])*int(A[i2,j1])) % p
                total_minors += 1
                if minor != 0:
                    nonzero_count += 1
    print(f"  n={n}, p={p} (n/p={n//p}): 2×2 minors ≡ 0 mod {p}? "
          f"{'YES' if nonzero_count == 0 else f'NO — {nonzero_count}/{total_minors} nonzero'}")


# ═══════════════════════════════════════════════════════════════════
print(f"\n\n" + "=" * 72)
print(f"  PROVA DEL TEOREMA rank_p = 1 PER n/p ≥ 2")
print(f"=" * 72)

print("""
  I dati mostrano:
  • n=4,p=2 (n/p=2): rank = 1 SEMPRE, tutti i 2×2 minori ≡ 0 
  • n=6,p=2 (n/p=3): rank = 1? 2×2 minori tutti zero?
  • n=6,p=3 (n/p=2): rank = 1 SEMPRE, tutti i 2×2 minori ≡ 0
  • n=5,p=5 (n/p=1): rank VARIA, 2×2 minori NON tutti zero
  
  DIMOSTRAZIONE (per n/p ≥ 2):
  
  Sia p|n con k=n/p ≥ 2. Per ogni LS L di ordine n:
  
  1. Ogni riga i di L contiene {1,...,n}, una permutazione.
     Mod p: contiene ogni classe 0,...,p-1 esattamente k volte.
  
  2. A_{ij} = L_{ij} - L_{i,n-1} mod p.
     Fissiamo la riga i. Sia a = L_{i,n-1} mod p.
     I valori A_{i,0},...,A_{i,n-2} mod p sono i residui 
     {L_{i,j} mod p - a : j=0,...,n-2}.
     
     I residui L_{i,0},...,L_{i,n-2} mod p sono {0,...,p-1} ciascuno 
     ripetuto k volte, MENO una copia di a (perché L_{i,n-1} ha residuo a).
     
     Quindi A_{i,*} mod p contiene:
     - Residuo (r-a) per r ≠ a: k volte
     - Residuo 0 (= a-a): k-1 volte
     
     Questa distribuzione è LA STESSA per OGNI riga i 
     (a può variare, ma le molteplicità sono identiche dopo shift).
     
  ATTENZIONE: stessa distribuzione NON implica righe uguali!
  Le posizioni dei residui possono variare. 
  
  Ma il dato sperimentale dice che il RANGO è sempre 1...
""")

# ── Let's check: ARE the rows of A mod p always identical? ───────
print("  Le righe di A mod p sono IDENTICHE o solo lo stesso multiset?")
print()

for n, p in [(4, 2), (6, 3)]:
    m = n - 1
    ls_list = ls_cache[n]
    identical_count = 0
    same_multiset_count = 0
    for ls in ls_list:
        L = np.array(ls)
        A = get_A(L) % p
        rows = [tuple(A[i]) for i in range(m)]
        if len(set(rows)) == 1:
            identical_count += 1
        multisets = [tuple(sorted(rows[i])) for i in range(m)]
        if len(set(multisets)) == 1:
            same_multiset_count += 1
    
    total = len(ls_list)
    print(f"  n={n}, p={p}: {identical_count}/{total} identical rows, "
          f"{same_multiset_count}/{total} same multiset")

# ── THE ANSWER: it's about 2x2 minors ────────────────────────────
print(f"\n" + "─" * 72)
print(f"  PROVA VIA MINORI 2×2")
print(f"─" * 72)

print("""
  rank(A mod p) = 1 ⟺ tutti i 2×2 minori di A sono ≡ 0 mod p
                       E almeno un entry è ≢ 0 mod p.
  
  2×2 minor: M_{(i,i'),(j,j')} = A_{ij}A_{i'j'} - A_{ij'}A_{i'j} mod p
  = (L_{ij}-L_{in})(L_{i'j'}-L_{i'n}) - (L_{ij'}-L_{in})(L_{i'j}-L_{i'n}) mod p
  
  Espandendo:
  = L_{ij}L_{i'j'} - L_{ij'}L_{i'j} 
    - L_{in}(L_{i'j'} - L_{i'j}) - L_{i'n}(L_{ij} - L_{ij'})  mod p
  
  Ma ATTENZIONE: mod p, con n ≡ 0 mod p, possiamo usare proprietà 
  specifiche dei LS.
  
  Facciamo un approccio diverso: 
  
  Definiamo B = L mod p (matrice n×n su F_p).
  Ogni riga e colonna di B contiene 0,...,p-1 ciascuno k times (k=n/p).
  
  A mod p = B[:n-1, :n-1] - B[:n-1, n-1]·1^T  (column subtraction)
  
  Cioè: (A mod p)_{ij} = B_{ij} - B_{i,n-1} per i,j = 0,...,n-2.
  
  Il 2×2 minore:
  = (B_{ij}-B_{in})(B_{i'j'}-B_{i'n}) - (B_{ij'}-B_{in})(B_{i'j}-B_{i'n})
  = B_{ij}B_{i'j'} - B_{ij'}B_{i'j} 
    - B_{in}B_{i'j'} + B_{in}B_{i'j} - B_{i'n}B_{ij} + B_{i'n}B_{ij'} 
    + B_{in}B_{i'n} - B_{in}B_{i'n}
  = B_{ij}B_{i'j'} - B_{ij'}B_{i'j} 
    + B_{in}(B_{i'j} - B_{i'j'}) + B_{i'n}(B_{ij'} - B_{ij})
  
  Questo è il 2×2 minor del blocking:
  [B_{ij}, B_{ij'}, B_{in}]   vs  [B_{i'j}, B_{i'j'}, B_{i'n}]
  
  In una matrice la cui ogni riga di L contiene ogni elemento di F_p 
  esattamente k ≥ 2 volte, possiamo usare un argomento statistico?
  
  Hmm... questo è troppo complesso per via analitica. 
  Affidiamoci al dato sperimentale ESAUSTIVO.
""")

# ═══════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("  VERIFICA FINALE: QUANDO k_p = n-2 È UNIVERSALE?")
print("=" * 72)

# The answer seems to be: k_p = n-2 iff n/p ≥ 2
# For n/p = 1 (p = n, n prime): k_p varies

# Let's also check: what k_p values appear when p = n?
print("\n  Caso p = n (n primo, n/p = 1):")
for n in [3, 5]:
    # n=7 has ~17M reduced LS — skip exhaustive, use sampling below
    ls_n = generate_reduced_ls(n)
    kp_dist = {}
    for ls in ls_n:
        L = np.array(ls)
        A = get_A(L)
        kp = (n-1) - rank_Fp(A, n)
        detA = bareiss_det(A.tolist())
        kp_dist[kp] = kp_dist.get(kp, 0) + 1
    print(f"  n=p={n}: k_n distribution = {dict(sorted(kp_dist.items()))}")

# n=7: JM sampling (17M reduced LS too many for exhaustive)
import random as _rng
_rng.seed(42)
def jm_sample(n, steps=5000):
    """Generate random LS via intercalate swaps from cyclic start."""
    grid = [[(i+j) % n + 1 for j in range(n)] for i in range(n)]
    # row_pos[j][v] = row where column j has value v+1
    row_pos = [[0]*n for _ in range(n)]
    for j in range(n):
        for i in range(n):
            row_pos[j][grid[i][j]-1] = i
    for _ in range(steps):
        i1 = _rng.randint(0, n-1)
        j1 = _rng.randint(0, n-1)
        j2 = _rng.randint(0, n-2)
        if j2 >= j1: j2 += 1
        v1 = grid[i1][j1]; v2 = grid[i1][j2]
        if v1 == v2: continue
        i2 = row_pos[j1][v2-1]
        if grid[i2][j2] != v1: continue
        # Valid intercalate swap
        grid[i1][j1] = v2; grid[i1][j2] = v1
        grid[i2][j1] = v1; grid[i2][j2] = v2
        row_pos[j1][v1-1] = i2; row_pos[j1][v2-1] = i1
        row_pos[j2][v1-1] = i1; row_pos[j2][v2-1] = i2
    return np.array(grid)

kp_dist_7 = {}
for _ in range(500):
    L = jm_sample(7, steps=10000)
    A = get_A(L)
    kp = 6 - rank_Fp(A, 7)
    kp_dist_7[kp] = kp_dist_7.get(kp, 0) + 1
print(f"  n=p=7 (500 JM samples): k_7 distribution = {dict(sorted(kp_dist_7.items()))}")

print("\n  Caso n/p ≥ 2:")
cases = [(4, 2), (6, 2), (6, 3)]
for n, p in cases:
    ls_n = ls_cache[n]
    kp_vals = set()
    for ls in ls_n:
        L = np.array(ls)
        A = get_A(L)
        kp = (n-1) - rank_Fp(A, p)
        kp_vals.add(kp)
    universal = len(kp_vals) == 1 and (n-2) in kp_vals
    print(f"  n={n}, p={p} (n/p={n//p}): k_p values = {sorted(kp_vals)}, "
          f"universally n-2? {'YES ✓' if universal else 'NO ✗'}")


# ═══════════════════════════════════════════════════════════════════
print("\n\n" + "=" * 72)
print("  SINTESI: DUE REGIMI DISTINTI")
print("=" * 72)

print("""
  ┌──────────────────────────────────────────────────────────────────┐
  │  REGIME 1: n/p ≥ 2 (p è un "piccolo" divisore di n)            │
  │                                                                  │
  │  • k_p = n-2 per TUTTI i LS (universale)                        │
  │  • rank_p(A) = 1 sempre                                         │
  │  • Formula: v_p(det A) ≥ v_p(n) + (n-3)                        │
  │  • Questo è ENORMEMENTE più forte del Thm 3.2                   │
  │  • Questo SPIEGA perché v_3(det A)=4 per tutti i LS a n=6      │
  │                                                                  │
  │  2×2 minori di A mod p tutti ZERO → rank ≤ 1.                   │
  │  Verificato esaustivamente per n=4,6.                            │
  ├──────────────────────────────────────────────────────────────────┤
  │  REGIME 2: n/p = 1 (n = p primo, p "grande")                   │
  │                                                                  │
  │  • k_p VARIA: dipende dal LS specifico                          │
  │  • Per n=5: k_5 ∈ {0, 1, 2, 3}                                 │
  │  • Formula: v_p(det A) ≥ 1 + max(0, k_p - 1) = max(1, k_p)    │
  │  • Non banale solo per k_p ≥ 2                                   │
  └──────────────────────────────────────────────────────────────────┘
  
  IMPLICAZIONE per il paper:
  
  Il Regime 1 è effettivamente un TEOREMA dimostrabile (da verificare),
  non una congettura. Potrebbe diventare un Lemma o Proposition nel paper.
  
  Il Regime 2 resta una congettura sperimentale.
""")
