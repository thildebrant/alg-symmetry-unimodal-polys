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

def make_latex_table(results):
    table_header = textwrap.dedent(r"""
    \begin{table}[h]
    \centering
    \begin{xltabular}{\textwidth}{@{}llp{4.5cm}p{4cm}@{}}
    \hline
    $(n,d)$ & \multicolumn{1}{c}{$I(C_n^d,x)$} & Modularity & Roots \\ \hline
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
        mod_str = r" \\ ".join(entry.get('mod_check', []))  # Double backslash for LaTeX
        mod_cell = r"\footnotesize$\begin{array}[t]{@{}l@{}}" + mod_str + r"\end{array}$"
        
        # Roots
        roots = [f"({r}, {mult})" for r, mult in entry.get('roots', [])]
        roots_cell = r"\footnotesize$\begin{array}[t]{@{}l@{}}" + r" \\ ".join(roots) + r"\end{array}$"
        
        row = rf"$({n},{d})$ & {poly_math} & {mod_cell} & {roots_cell} \\\\ \hline"
        table_rows.append(row)
    
    table_footer = textwrap.dedent(r"""
    \hline
    \end{xltabular}
    \caption{Summary of computations for selected $(n,d)$}
    \end{table}
    """).strip()
    
    return "\n".join([table_header] + table_rows + [table_footer])

if __name__ == "__main__":
    nd_list = [(3,1), (3,2), (3,3), (4,1), (4,2), (4,3), (5,1), (5,2), (5,3), (6,1), (6,2), (6,3), (7,1), (7,2), (7,3), (8,1), (8,2),(9,1),(9,2)]
    results = analyze_polynomials(nd_list, compute_roots=True, compute_galois=False, check_modular=True)
    latex_code = make_latex_table(results)
    print(latex_code)
