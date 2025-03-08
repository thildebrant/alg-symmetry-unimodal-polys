from sage.all import *

def analyze_Cn_strong_product(n, d, p=None):
    r"""
    Analyze the independence polynomial I(C_n^d, x) in Sage, covering:
    
      (1) Modular collapse checks:
          - If n is prime, verify I(x) ≡ 1 (mod n).
          - If n is composite, for each prime p | n (or a user-specified p),
            verify I(x) ≡ (x+1)^(p^d) mod p.
            
      (2) Factorization and splitting field:
          - Compute I(x) in QQ[x], factor over QQ,
            then factor over its splitting field K.
            
      (3) Galois group or automorphism group:
          - If I(x) is irreducible in QQ[x], call galois_group().
          - If reducible, compute Galois groups of irreducible factors (degree>1).
          - Attempt to get automorphism group of K if possible.
            
    Note: This function can be expensive for larger n, d, or polynomials above degree ~10.
    """

    # Construct the d-fold strong product of C_n.
    G = graphs.CycleGraph(n)
    for _ in range(d-1):
        G = G.strong_product(graphs.CycleGraph(n))

    # Complement G => clique polynomial => independence polynomial of G.
    Gc = G.complement()
    I_poly = Gc.clique_polynomial().change_variable_name('x')
    
    print(f"\n=== Analyzing C_{n}^{d} ===")
    print(f"Independence polynomial I(x): {I_poly}")
    
    # Check modular collapse.
    if n.is_prime():
        # Check I(x) mod n == 1
        I_mod_n = I_poly.change_ring(GF(n))
        is_one = (I_mod_n == I_mod_n.parent()(1))
        print(f"\n(1) Modular check (prime n={n}):  I(x) ≡ 1 mod {n}? -> {is_one}")
        print(f"    I(x) mod {n} = {I_mod_n}")
    else:
        # Composite n => check each prime divisor (or use user-specified p).
        pdvs = n.prime_divisors()
        if not pdvs:
            print(f"\nWarning: n={n} has no prime divisors?? (Should not happen unless n=1).")
            return
        
        print(f"\n(1) Modular check (composite n={n}):")
        if p is not None:
            # Just check the user-specified prime.
            pdvs = [p] if p in pdvs else []
            if not pdvs:
                print(f"  Error: specified p={p} is not a divisor of n={n}.")
                return
        
        for prime_div in pdvs:
            I_mod_p = I_poly.change_ring(GF(prime_div))
            Rp = GF(prime_div)['X']
            X = Rp.gen()
            # Use prime_div**d
            target = (X + 1)^(prime_div**d)

            eq_check = (Rp(I_mod_p) == target)
            print(f"  - p={prime_div}:  I(x) ≡ (x+1)^{prime_div**d} mod {prime_div}? -> {eq_check}")
            if not eq_check:
                print(f"    => I(x) mod {prime_div} is {Rp(I_mod_p)}")

    
    # Factorization in QQ.
    R.<x> = QQ[]
    I_QQ = R(I_poly)  # interpret polynomial in QQ[x].
    
    print("\n(2) Factorization over QQ:")
    fact_qq = I_QQ.factor()
    print(f"    {fact_qq}")
    
    # Splitting field and factor over that field.
    print("\n(3) Splitting field computation:")
    try:
        K = I_QQ.splitting_field('alpha')
        print(f"    K = {K}")
        Ifact = I_QQ.change_ring(K).factor()
        print(f"    Factorization in K:\n      {Ifact}")
    except (RuntimeError, ValueError, NotImplementedError) as e:
        print(f"    Splitting field computation not feasible: {e}")
        return
    
    # Galois group or automorphism group analysis.
    print("\n(4) Galois/Automorphism group analysis:")
    if I_QQ.is_irreducible():
        try:
            Galois_grp = I_QQ.galois_group()
            grp_str = Galois_grp.structure_description()
            print(f"    Galois group (irreducible case) = {grp_str}")
        except (NotImplementedError, ValueError, RuntimeError) as e:
            print(f"    Galois group computation failed: {e}")
    else:
        print("    Polynomial is reducible over QQ.")
        factors = I_QQ.factor()
        for fct, expt in factors:
            if fct.degree() > 1:
                try:
                    grp = fct.galois_group()
                    print(f"      Factor: {fct}, deg={fct.degree()}, Galois group => {grp.structure_description()}")
                except (NotImplementedError, ValueError, RuntimeError) as e:
                    print(f"      Galois group for {fct} not feasible: {e}")
        # Attempt full splitting field automorphisms:
        try:
            auto_grp = K.automorphism_group()
            print(f"    Automorphism group of K:\n      {auto_grp}")
        except AttributeError:
            auts = K.automorphisms()
            print(f"    List of automorphisms in K:\n      {auts}")


# Example usage:
if __name__ == "__main__":
    # Prime case: n=5, d=1.
    analyze_Cn_strong_product(5, 1)
    
    # Composite case: n=6, d=1.
    analyze_Cn_strong_product(6, 1)

# Example output
'''
\begin{lstlisting}
=== Analyzing C_5^1 ===
Independence polynomial I(x): 5*x^2 + 5*x + 1

(1) Modular check (prime n=5):  I(x) ≡ 1 mod 5? -> True
    I(x) mod 5 = 1

(2) Factorization over QQ:
    (5) * (x^2 + x + 1/5)

(3) Splitting field computation:
    K = Number Field in alpha with defining polynomial x^2 + 5*x + 5
    Factorization in K:
      (5) * (x - 1/5*alpha) * (x + 1/5*alpha + 1)

(4) Galois/Automorphism group analysis:
    Galois group (irreducible case) = C2

=== Analyzing C_6^1 ===
Independence polynomial I(x): 2*x^3 + 9*x^2 + 6*x + 1

(1) Modular check (composite n=6):
  - p=2:  I(x) ≡ (x+1)^2 mod 2? -> True
  - p=3:  I(x) ≡ (x+1)^3 mod 3? -> False
    => I(x) mod 3 is 2*X^3 + 1

(2) Factorization over QQ:
    (2) * (x + 1/2) * (x^2 + 4*x + 1)

(3) Splitting field computation:
    K = Number Field in alpha with defining polynomial x^2 + 4*x + 1
    Factorization in K:
      (2) * (x - alpha) * (x + 1/2) * (x + alpha + 4)

(4) Galois/Automorphism group analysis:
    Polynomial is reducible over QQ.
      Factor: x^2 + 4*x + 1, deg=2, Galois group => C2
    List of automorphisms in K:
      [
Ring endomorphism of Number Field in alpha with defining polynomial x^2 + 4*x + 1
  Defn: alpha |--> alpha,
Ring endomorphism of Number Field in alpha with defining polynomial x^2 + 4*x + 1
  Defn: alpha |--> -alpha - 4
]
\end{lstlisting}
'''
