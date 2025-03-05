import logging
from functools import lru_cache
import sys
import networkx as nx
import numpy as np
import psutil

memory_info = psutil.virtual_memory()
print(f"Total memory: {memory_info.total / (1024 ** 3):.2f} GB")
print(f"Available memory: {memory_info.available / (1024 ** 3):.2f} GB")


# Logging setup
logging.basicConfig(
    filename="debug_log_recursive.txt",
    filemode="w",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

# Configure recursion limit (if needed)
sys.setrecursionlimit(3000)
logging.info(f"Recursion limit set to {sys.getrecursionlimit()}")

# Counter for the number of recursions
recursion_count = 0

def precompute_neighbors(graph):
    """Precompute and cache neighbors for each node."""
    return {node: set(graph.neighbors(node)) for node in graph.nodes()}

def independence_polynomial(graph):
    global recursion_count
    recursion_count = 0  # Reset the counter

    neighbors_dict = precompute_neighbors(graph)  

    @lru_cache(maxsize=50000000)  # Increase cache size
    def _independent_polynomial_recursive(nodes):
        global recursion_count
        recursion_count += 1
        if recursion_count % 10000000 == 0:
            logging.info(f"Recursions: {recursion_count}")

        if not nodes:
            return np.array([1], dtype=np.int64)
        
        node = nodes[0]
        remaining_nodes = nodes[1:]

        exclude_poly = _independent_polynomial_recursive(tuple(remaining_nodes))

        neighbors = neighbors_dict[node]
        remaining_after_neighbors = tuple(n for n in remaining_nodes if n not in neighbors)

        include_poly = _independent_polynomial_recursive(remaining_after_neighbors)

        result_poly_size = max(len(exclude_poly), len(include_poly) + 1)
        result_poly = np.zeros(result_poly_size, dtype=np.int64)
        result_poly[:len(exclude_poly)] += exclude_poly
        result_poly[1:len(include_poly)+1] += include_poly

        return result_poly

    nodes = tuple(graph.nodes())
    coefficients = _independent_polynomial_recursive(nodes)

    poly_str = " + ".join(
        f"{coeff}x^{i}" if i > 1 else (f"{coeff}x" if i == 1 else f"{coeff}")
        for i, coeff in enumerate(coefficients)
        if coeff > 0
    )
    logging.info(f"Final polynomial: {poly_str}")
    logging.info(f"Total recursions: {recursion_count}")
    logging.info(f"Cache info: {_independent_polynomial_recursive.cache_info()}")
    return poly_str

def compute_graph_recursive(n, k):
    g = nx.cycle_graph(n)
    for _ in range(1, k):
        g = nx.strong_product(g, nx.cycle_graph(n))

    poly = independence_polynomial(g)
    return poly

if __name__ == "__main__":
    n = 5
    k = 3
    poly = compute_graph_recursive(n, k)
    print(f"Independence Polynomial: {poly}")
