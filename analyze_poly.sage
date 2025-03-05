############################################################
# File: analyze_poly.sage
#
# Sage script to:
#   1) Build or read independence polynomials I(C_n^{d}, x).
#   2) Prime-factorize coefficients.
#   3) Perform modular checks (if desired).
#   4) Compute Galois group or numerical roots (optional).
#   5) Generate LaTeX table summarizing data.
#
# Usage inside Sage:
#    sage: load("analyze_poly.py")
#    sage: results = analyze_polynomials([(3,1),(4,2),...])  # or from your stored polynomials
#    sage: print(make_latex_table(results))
############################################################

from sage.all import *
import textwrap

def analyze_polynomials(
    nd_list,
    from_strings=None,
    check_modular=True,
    compute_galois=False,
    compute_roots=False
):
    r"""
    Analyze independence polynomials for each (n,d) in nd_list.
    
    If from_strings is not None, it should be a dict or list
    providing polynomial strings keyed by (n,d).
    Otherwise, we construct I(C_n^d, x) from Sage's graph routines.
    
    Returns a list of dictionaries, each storing:
       {
         'n': n,
         'd': d,
         'poly_s':  string form of the polynomial,
         'coeff_factors': {k: {prime:exponent, ...}, ...}   # prime factorization of each coeff
         'mod_check': [ details of modular checks, if check_modular=True ],
         'galois_group': "...",
         'roots': [approx roots],   # if compute_roots
         ...
       }
    """
    results = []
    
    for (n,d) in nd_list:
        data = {'n': n, 'd': d}
        
        # (A) Build or read the polynomial
        if from_strings and (n,d) in from_strings:
            poly_s = from_strings[(n,d)]
            R.<x> = ZZ[]
            try:
                poly = R(poly_s)
            except:
                print(f"Could not parse polynomial string for (n={n}, d={d}).")
                continue
        else:
            # Construct polynomial from Sage Graph approach
            G = graphs.CycleGraph(n)
            for _ in range(d-1):
                G = G.strong_product(graphs.CycleGraph(n))
            # independence polynomial = clique polynomial of complement
            Gc = G.complement()
            poly_sage = Gc.clique_polynomial().change_variable_name('x')
            # convert to a polynomial over ZZ
            R.<x> = ZZ[]
            poly = R(poly_sage)
            poly_s = str(poly)
        
        data['poly_s'] = poly_s
        
        # (B) Prime factorization of coefficients
        coeff_factor_dict = {}
        coeffs = poly.coefficients(sparse=False)
        # Sage returns coefficients in descending order if sparse=False,
        # but let's confirm index => power of x
        # Actually, polynomial.degree() = len(coeffs)-1 for usual representation
        deg = poly.degree()
        # coeff[i] is the coefficient of x^(deg - i) if we use standard descending indexing
        # We'll unify by reversing or just re-index carefully:
        for i, c in enumerate(reversed(coeffs)):
            # i = 0 => x^0, i = 1 => x^1, etc.
            if c != 0:
                factor_dict = dict(factor(c))
                coeff_factor_dict[i] = factor_dict
        data['coeff_factors'] = coeff_factor_dict
        
        # (C) Optional: Modular checks for prime/composite n
        if check_modular:
            mod_details = []
            if n.is_prime():
                # check I(x) mod n == 1
                pmod = poly.change_ring(GF(n))
                is_one = bool(pmod == pmod.parent()(1))
                mod_details.append(f"I(x) ≡ 1 (mod {n})? => {is_one}")
            else:
                # Check for each prime divisor
                pdvs = n.prime_divisors()
                for p in pdvs:
                    pmod = poly.change_ring(GF(p))
                    # We want to see if pmod == (x+1)^(p^d)
                    Rp.<X> = GF(p)[]
                    target = (X+1)^(p**d)
                    eq_flag = (Rp(pmod) == target)
                    mod_details.append(f"I(x) ≡ (x+1)^{p**d} mod {p}? => {eq_flag}")
            data['mod_check'] = mod_details
        
        # (D) Optional: Galois group
        if compute_galois:
            poly_QQ = poly.change_ring(QQ)
            if poly_QQ.degree() <= 1:
                data['galois_group'] = "Trivial or linear"
            else:
                try:
                    if poly_QQ.is_irreducible():
                        Ggrp = poly_QQ.galois_group()
                        data['galois_group'] = Ggrp.structure_description()
                    else:
                        # multiple factors => do a partial analysis
                        fct_str = []
                        for fct,ex in poly_QQ.factor():
                            if fct.degree() > 1:
                                try:
                                    gsub = fct.galois_group().structure_description()
                                    fct_str.append(f"{fct} => {gsub}")
                                except Exception as e:
                                    fct_str.append(f"{fct} => Galois group fail {e}")
                        data['galois_group'] = "; ".join(fct_str) if fct_str else "All linear factors"
                except Exception as e:
                    data['galois_group'] = f"Galois group fail: {e}"
        
        # (E) Optional: numeric roots
        if compute_roots:
            # approximate numeric roots (over RR or CC)
            poly_RR = poly.change_ring(RR)
            # We can find approximate roots with .roots() in RR, but if some are complex,
            # we do .roots() in CC
            try:
                rts = poly_RR.roots(ring=CC)
                # each item is (root_value, multiplicity)
                data['roots'] = [(str(r[0]), r[1]) for r in rts]
            except Exception as e:
                data['roots'] = [f"Root finding fail: {e}"]
        
        results.append(data)
    
    return results


def make_latex_table(results):
    r"""
    Generate a LaTeX tabular environment summarizing the analysis.
    
    We produce columns for (n, d), polynomial, prime factorization of coefficients,
    possibly Galois group info, etc.
    """
    # You can adapt columns as you wish. Let's do something basic:
    #   n & d & "Polynomial" & "Coefficient prime factors" & "Modular check" & "Galois group"
    
    table_header = textwrap.dedent(r"""
    \begin{table}[h]
    \centering
    \begin{tabular}{lllll}
    \hline
    $(n,d)$ & $I(C_n^d,x)$ & Coeff. Factorization & Mod. Check & Galois \\ \hline
    """).strip()
    
    table_rows = []
    
    for entry in results:
        n, d = entry['n'], entry['d']
        poly_s = entry.get('poly_s', '')
        
        # Format coefficient factorization
        # We'll just show partial info or you can do a more elaborate approach
        cfac = entry.get('coeff_factors', {})
        # cfac is a dict: power -> {prime: exp, ...}
        # We'll format it as e.g. x^0:2^3 * 5  ; x^1: ...
        cfac_parts = []
        for k in sorted(cfac.keys()):
            prime_exp = []
            for p, e in sorted(cfac[k].items()):
                if e == 1:
                    prime_exp.append(str(p))
                else:
                    prime_exp.append(f"{p}^{e}")
            cfac_parts.append(f"x^{k}: " + "*".join(prime_exp))
        cfac_str = r"\\".join(cfac_parts)  # linebreak in latex cell
        
        # Modular check
        modcheck = entry.get('mod_check', [])
        modcheck_str = r"\\".join(modcheck)
        
        # Galois group
        gal = entry.get('galois_group',"")
        
        row = f"{(n,d)} & {poly_s} & \\tiny {cfac_str} & \\tiny {modcheck_str} & \\tiny {gal} \\\\"
        table_rows.append(row)
    
    table_footer = textwrap.dedent(r"""
    \hline
    \end{tabular}
    \caption{Summary of computations for selected $(n,d)$}
    \end{table}
    """).strip()
    
    latex_table = "\n".join([table_header] + table_rows + [table_footer])
    return latex_table

#####################################
# Example usage if running in Sage:
#####################################
if __name__ == "__main__":
    # We can define (n,d) pairs
    nd_list = [(3,1), (4,2), (5,2), (6,2)]
    
    # or we can supply a dictionary of polynomials you already have
    # from_strings = {
    #     (3,1): "3*x+1",
    #     (4,2): "12*x^4+48*x^3+56*x^2+16*x+1",
    #     ...
    # }
    # results = analyze_polynomials(nd_list, from_strings=from_strings, compute_roots=True, compute_galois=True)
    
    # For demonstration, let's build them from Sage's cycle-product approach
    results = analyze_polynomials(nd_list, from_strings=None, compute_roots=True, compute_galois=True)
    
    # Produce LaTeX table
    latex_code = make_latex_table(results)
    print(latex_code)
