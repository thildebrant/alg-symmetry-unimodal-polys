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

def analyze_polynomials(nd_list, from_strings=None, check_modular=True,
                       compute_galois=True, compute_roots=True):
    results = []
    
    for (n,d) in nd_list:
        data = {'n': n, 'd': d}
        poly = None  # Initialize poly
        
        # Polynomial generation/parsing
        if from_strings and (n,d) in from_strings:
            try:
                R.<x> = ZZ[]
                poly = R(from_strings[(n,d)])
            except Exception as e:
                print(f"Error parsing ({n},{d}): {str(e)}")
                continue
        else:
            try:
                G = graphs.CycleGraph(n)
                for _ in range(d-1):
                    G = G.strong_product(graphs.CycleGraph(n))
                Gc = G.complement()
                poly = Gc.clique_polynomial().change_variable_name('x')
            except Exception as e:
                print(f"Error generating ({n},{d}): {str(e)}")
                continue
        
        if poly is None:  # Safety check
            continue
            
        # Format polynomial
        poly_s = latex(poly)
        terms = poly_s.split(' + ')
        poly_s = ' + \\\\\n'.join([' + '.join(terms[i:i+3]) for i in range(0, len(terms), 3)])
        data['poly_s'] = f"\\begin{{aligned}}{poly_s}\\end{{aligned}}"
        
        # Modular checks
        if check_modular:
            mod_details = []
            if n.is_prime():
                pmod = poly.change_ring(GF(n))
                is_one = bool(pmod == pmod.parent()(1))
                mod_details.append(f"I(x) \\equiv 1 \\mod {n}? \\Rightarrow {is_one}")
            else:
                for p in n.prime_divisors():
                    pmod = poly.change_ring(GF(p))
                    target = (poly.parent().gen() + 1)**(p**d)
                    eq_flag = (pmod == target)
                    mod_details.append(f"I(x) \\equiv (x+1)^{{{p**d}}} \\mod {p}? \\Rightarrow {eq_flag}")
            data['mod_check'] = mod_details
        
        # Galois group
        if compute_galois:
            try:
                poly_QQ = poly.change_ring(QQ)
                if poly_QQ.degree() <= 1:
                    data['galois_group'] = "Trivial"
                else:
                    if poly_QQ.is_irreducible():
                        G = poly_QQ.galois_group()
                        data['galois_group'] = G.structure_description()
                    else:
                        factors = []
                        for f, _ in poly_QQ.factor():
                            if f.degree() > 1:
                                try:
                                    G = f.galois_group()
                                    factors.append(f"{latex(f)} \\Rightarrow {G.structure_description()}")
                                except:
                                    factors.append(f"{latex(f)} \\Rightarrow Error")
                        data['galois_group'] = "; ".join(factors)
            except Exception as e:
                data['galois_group'] = f"Error: {str(e)}"
        
        # Roots
        if compute_roots:
            try:
                poly_CC = poly.change_ring(CC)
                roots = poly_CC.roots(multiplicities=True)
                formatted_roots = []
                for r, mult in roots:
                    if abs(r.imag()) < 1e-6:
                        formatted_roots.append((f"{r.real():.5f}", mult))
                    else:
                        formatted_roots.append((f"{r.real():.5f} \\pm {abs(r.imag()):.5f}i", mult))
                data['roots'] = formatted_roots
            except Exception as e:
                data['roots'] = [("Error", 1)]
        
        results.append(data)
    
    return results
# Keep the make_latex_table function from previous working version
def make_latex_table(results):
    table_header = textwrap.dedent(r"""
    \begin{table}[h]
    \centering
    \begin{xltabular}{\textwidth}{@{}llXcc@{}}
    \hline
    $(n,d)$ & \multicolumn{1}{c}{$I(C_n^d,x)$} & Mod. Check & Galois & Roots \\ \hline
    """).strip()
    
    table_rows = []
    
    for entry in results:
        if not entry:  # Skip empty entries
            continue
            
        n, d = entry['n'], entry['d']
        poly_s = entry.get('poly_s', '')
        
        # Polynomial formatting
        poly_math = r"\footnotesize$\begin{aligned}" + poly_s + r"\end{aligned}$"
        
        # Modular checks
        mod_str = r" \ ".join(entry.get('mod_check', []))  # Double backslash for LaTeX
        mod_cell = r"\footnotesize$\begin{array}[t]{@{}l@{}}" + mod_str + r"\end{array}$"
        
        # Galois group
        gal_str = r" \ ".join(entry.get('galois_group', "").split("; "))
        gal_cell = r"\footnotesize$\begin{array}[t]{@{}l@{}}" + gal_str + r"\end{array}$"
        
        # Roots
        roots = [f"({r}, {mult})" for r, mult in entry.get('roots', [])]
        roots_cell = r"\footnotesize$\begin{array}[t]{@{}l@{}}" + r" \\ ".join(roots) + r"\end{array}$"
        
        row = f"$({n},{d})$ & {poly_math} & {mod_cell} & {gal_cell} & {roots_cell} \\\\"
        table_rows.append(row)
    
    table_footer = textwrap.dedent(r"""
    \hline
    \end{xltabular}
    \caption{Summary of computations for selected $(n,d)$}
    \end{table}
    """).strip()
    
    return "\n".join([table_header] + table_rows + [table_footer])
def format_modcheck(mod_details):
    if not mod_details:
        return ""
    lines = r" \\ ".join([m.replace("=>", r"\Rightarrow") for m in mod_details])
    return r"\footnotesize$\begin{array}[t]{@{}l@{}} " + lines + r" \end{array}$"

def format_galois(gal_data):
    if ";" in gal_data:
        parts = gal_data.split("; ")
        lines = r" \\ ".join(parts)
        return r"\footnotesize$\begin{array}[t]{@{}l@{}} " + lines + r" \end{array}$"
    return r"\footnotesize$" + gal_data + r"$"

def format_roots(roots_list):
    if not roots_list:
        return ""
    lines = r" \\ ".join([f"({r}, {mult})" for r, mult in roots_list])
    return r"\footnotesize$\begin{array}[t]{@{}l@{}} " + lines + r" \end{array}$"

if __name__ == "__main__":
    nd_list = [(3,1), (4,2), (5,2), (6,2)]
    results = analyze_polynomials(nd_list, compute_roots=True, compute_galois=True, check_modular=True)    
    latex_code = make_latex_table(results)
    print(latex_code)
