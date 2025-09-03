import numpy as np
import networkx as nx
from collections import defaultdict
import math

def solve_matrix_decomposition(A, k):
    """
    Solves the matrix decomposition problem using a graph-based approach.

    Args:
        A (np.ndarray): The input binary matrix (m x n).
        k (int): The number of initial rows with placement constraints.

    Returns:
        dict: A dictionary containing the results of the optimization.
    """
    m, n = A.shape
    print(f"Solving for a {m}x{n} matrix with k={k}\n")

    if m == 0 or n == 0:
        return {
            "optimal_cost": 0,
            "U_rows_count": 0,
            "L_max_rows": 0,
            "U_rows": set(),
            "U_cols": set(),
            "L_blocks": []
        }

    # 1. Build the bipartite graph
    G = nx.Graph()
    rows = [f"r{i}" for i in range(m)]
    cols = [f"c{i}" for i in range(n)]
    G.add_nodes_from(rows, bipartite=0)
    G.add_nodes_from(cols, bipartite=1)

    non_zero_indices = np.argwhere(A == 1)
    for r_idx, c_idx in non_zero_indices:
        G.add_edge(f"r{r_idx}", f"c{c_idx}")

    # 2. Find connected components
    connected_components = list(nx.connected_components(G))

    # 3. Process components to get row/column sets and sizes
    components_info = []
    for i, comp_nodes in enumerate(connected_components):
        comp_rows = {int(node[1:]) for node in comp_nodes if node.startswith('r')}
        comp_cols = {int(node[1:]) for node in comp_nodes if node.startswith('c')}
        if comp_rows: # Only consider components with rows
            components_info.append({
                "id": i,
                "rows": comp_rows,
                "cols": comp_cols,
                "row_count": len(comp_rows),
            })

    # 4. Separate special and regular components
    special_row_indices = set(range(k))
    special_component_rows = set()
    special_component_cols = set()
    regular_components = []

    # Find all components touching the first k rows and merge them
    special_comp_ids_to_merge = set()
    for comp in components_info:
        if not comp["rows"].isdisjoint(special_row_indices):
            special_comp_ids_to_merge.add(comp["id"])
            special_component_rows.update(comp["rows"])
            special_component_cols.update(comp["cols"])

    special_component = None
    if special_component_rows:
        special_component = {
            "rows": special_component_rows,
            "cols": special_component_cols,
            "row_count": len(special_component_rows)
        }

    for comp in components_info:
        if comp["id"] not in special_comp_ids_to_merge:
            regular_components.append(comp)

    # Sort regular components by row count, descending (for packing heuristic)
    regular_components.sort(key=lambda x: x["row_count"], reverse=True)

    def get_packing_cost(target_max_l_rows, special_comp, regular_comps):
        """
        Calculates the total cost (U_rows + max_L_rows) for a given target.
        Uses a First-Fit-Decreasing heuristic for bin packing.
        """
        
        # --- Case 1: Special component is forced into U ---
        u_rows_case1 = set()
        l_blocks_case1 = []
        
        if special_comp:
            u_rows_case1.update(special_comp["rows"])

        # Pack regular components
        bins_case1 = [] # Stores row counts of L blocks
        for comp in regular_comps:
            if comp["row_count"] > target_max_l_rows:
                u_rows_case1.update(comp["rows"])
                continue
            
            # Find first bin that fits
            placed = False
            for i in range(len(bins_case1)):
                if bins_case1[i] + comp["row_count"] <= target_max_l_rows:
                    bins_case1[i] += comp["row_count"]
                    placed = True
                    break
            if not placed:
                bins_case1.append(comp["row_count"])
        
        cost1 = len(u_rows_case1) + target_max_l_rows

        # --- Case 2: Special component is placed in L (if possible) ---
        cost2 = float('inf')
        if special_comp and special_comp["row_count"] <= target_max_l_rows:
            u_rows_case2 = set()
            l_blocks_case2 = []
            
            bins_case2 = [special_comp["row_count"]] # Start L with special comp

            # Pack regular components
            for comp in regular_comps:
                if comp["row_count"] > target_max_l_rows:
                    u_rows_case2.update(comp["rows"])
                    continue
                
                placed = False
                for i in range(len(bins_case2)):
                    if bins_case2[i] + comp["row_count"] <= target_max_l_rows:
                        bins_case2[i] += comp["row_count"]
                        placed = True
                        break
                if not placed:
                    bins_case2.append(comp["row_count"])
            
            cost2 = len(u_rows_case2) + target_max_l_rows
        
        return min(cost1, cost2)

    # 5. Binary search for the optimal max_L_rows (T)
    low = 0
    high = m  # Max possible value for max_L_rows is all rows
    optimal_t = m
    min_cost = m

    while low <= high:
        t = (low + high) // 2
        if t < 0: break

        cost_t = get_packing_cost(t, special_component, regular_components)
        
        # Check if a smaller T could be even better
        cost_t_minus_1 = get_packing_cost(t - 1, special_component, regular_components) if t > 0 else float('inf')

        if cost_t <= min_cost:
            min_cost = cost_t
            optimal_t = t

        if cost_t <= cost_t_minus_1:
            # It's getting cheaper or staying same as T increases, so maybe we can do better with an even smaller T
            high = t - 1
        else:
            # Cost increased when we went from T-1 to T, so optimal must be >= T
            low = t + 1

    # 6. Reconstruct the final optimal assignment with optimal_t
    
    # Case 1 final reconstruction
    u_rows_final1 = set()
    l_blocks_final1 = []
    if special_component:
        u_rows_final1.update(special_component["rows"])
    
    bins1 = []
    for comp in regular_components:
        if comp["row_count"] > optimal_t:
            u_rows_final1.update(comp["rows"])
        else:
            placed = False
            for i in range(len(bins1)):
                if bins1[i]["row_count"] + comp["row_count"] <= optimal_t:
                    bins1[i]["components"].append(comp)
                    bins1[i]["row_count"] += comp["row_count"]
                    placed = True
                    break
            if not placed:
                bins1.append({"components": [comp], "row_count": comp["row_count"]})
    cost_final1 = len(u_rows_final1) + optimal_t

    # Case 2 final reconstruction
    cost_final2 = float('inf')
    u_rows_final2 = set()
    l_blocks_final2 = []
    if special_component and special_component["row_count"] <= optimal_t:
        bins2 = [{"components": [special_component], "row_count": special_component["row_count"]}]
        for comp in regular_components:
            if comp["row_count"] > optimal_t:
                u_rows_final2.update(comp["rows"])
            else:
                placed = False
                for i in range(len(bins2)):
                    if bins2[i]["row_count"] + comp["row_count"] <= optimal_t:
                        bins2[i]["components"].append(comp)
                        bins2[i]["row_count"] += comp["row_count"]
                        placed = True
                        break
                if not placed:
                    bins2.append({"components": [comp], "row_count": comp["row_count"]})
        cost_final2 = len(u_rows_final2) + optimal_t
    
    # Choose the best final configuration
    if cost_final1 <= cost_final2:
        final_u_rows = u_rows_final1
        final_l_bins = bins1
    else:
        final_u_rows = u_rows_final2
        final_l_bins = bins2

    # Assemble final L blocks from bins
    final_l_blocks = []
    for i, bin_info in enumerate(final_l_bins):
        block_rows = set()
        block_cols = set()
        for comp in bin_info["components"]:
            block_rows.update(comp["rows"])
            block_cols.update(comp["cols"])
        final_l_blocks.append({
            "block_id": i,
            "row_indices": block_rows,
            "col_indices": block_cols,
            "row_count": len(block_rows)
        })
    
    all_l_rows = set().union(*(b['row_indices'] for b in final_l_blocks))
    all_l_cols = set().union(*(b['col_indices'] for b in final_l_blocks))
    
    # U contains all rows/cols not in L
    final_u_rows.update(set(range(m)) - all_l_rows)
    final_u_cols = set(range(n)) - all_l_cols

    return {
        "optimal_cost": len(final_u_rows) + optimal_t,
        "U_rows_count": len(final_u_rows),
        "L_max_rows": optimal_t,
        "U_rows": sorted(list(final_u_rows)),
        "U_cols": sorted(list(final_u_cols)),
        "L_blocks": final_l_blocks
    }


def main():
    """Main function to run a test case."""
    # --- Test Case Generation ---
    # You can change these parameters
    m_test = 13  # Number of rows
    n_test = 49  # Number of columns
    k_test = 1   # First k rows are special
    
    tri49 = np.array([
    # even-weight rows (rows 1 – 13) – this is G₀ in App. B of Bravyi-Haah
    list(map(int, "1111111111111110101010101010101010101010101010101")),
    list(map(int, "0000000000000000000111100110011000011001100110011")),
    list(map(int, "0000000000000001100000011001100110000000000000000")),
    list(map(int, "0000000000000000000000000000000001111000000001111")),
    list(map(int, "0000000000000000011110000000000000000111100000000")),
    list(map(int, "0000000000000000000001111000011110000000000000000")),
    list(map(int, "0000000000000000000000000111111110000000000000000")),
    list(map(int, "0000000000000000000000000000000001111111100000000")),
    list(map(int, "0000000000000000000000000000000000000000011111111")),
    list(map(int, "1010101010101010000000000000000000000000000000000")),
    list(map(int, "0110011001100110000000000000000000000000000000000")),
    list(map(int, "0001111000011110000000000000000000000000000000000")),
    list(map(int, "0000000111111110000000000000000000000000000000000"))
    ], dtype=int)

    hex_values = [
    "000000000000ffff",
    "000f000f000f000f",
    "1111111111111111",
    "000000ff000000ff",
    "0303030303030303",
    "0000000055555555",
    "ffffffffffffffff",
    "00000000ffffffff",
    "0000ffff0000ffff",
    "00ff00ff00ff00ff",
    "0f0f0f0f0f0f0f0f",
    "3333333333333333",
    "5555555555555555",
    "0000000000ff00ff",
    "0033003300330033",
    "0000000033333333",
    "0000000000330033",
    ]


    a64to2ccz = np.array([
    list(map(int, bin(int(h, 16))[2:].zfill(64)))
    for h in hex_values
    ], dtype=np.uint8)
 
    A_test = a64to2ccz


    print("--- Input Matrix (A) ---")
    # print(A_test)
    print(f"Shape: {A_test.shape}")
    print(f"Number of '1's: {np.sum(A_test)}")
    print("-" * 26)

    # --- Solve the problem ---
    result = solve_matrix_decomposition(A_test, k_test)

    # --- Print Results ---
    print("\n--- Optimization Results ---")
    print(f"Optimal Cost (U_rows + max_L_rows): {result['optimal_cost']}")
    print(f"  - Rows in U: {result['U_rows_count']}")
    print(f"  - Max rows in an L block: {result['L_max_rows']}")
    print("-" * 28)

    print(f"\nUpper Partition (U):")
    print(f"  - Row indices: {result['U_rows']}")
    print(f"  - Column indices: {result['U_cols']}")

    print(f"\nLower Partition (L) has {len(result['L_blocks'])} blocks:")
    for block in result['L_blocks']:
        print(f"  - L Block {block['block_id']}:")
        print(f"    - Row count: {block['row_count']}")
        print(f"    - Row indices: {sorted(list(block['row_indices']))}")
        print(f"    - Col indices: {sorted(list(block['col_indices']))}")

if __name__ == "__main__":
    main()

