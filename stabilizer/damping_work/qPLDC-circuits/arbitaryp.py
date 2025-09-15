import numpy as np
from itertools import permutations, combinations
import copy

class MatrixPartitionSolver:
    def __init__(self, matrix, k):
        """
        Initialize the solver with a binary matrix and constraint parameter k.
        
        Args:
            matrix: 2D numpy array with binary values (0 or 1)
            k: number of first rows that can only go to U or first row partition of L
        """
        self.original_matrix = np.array(matrix)
        self.k = k
        self.m, self.n = self.original_matrix.shape
        self.best_solution = None
        self.best_objective = float('inf')
        
    def check_anti_diagonal_structure(self, matrix_L, p, row_partitions, col_partitions):
        """
        Check if the lower half L satisfies the anti-diagonal block structure.
        
        Args:
            matrix_L: the lower half matrix
            p: number of partitions
            row_partitions: list of row partition sizes
            col_partitions: list of column partition sizes
            
        Returns:
            bool: True if structure is satisfied
        """
        # Get block boundaries
        row_bounds = [0] + [sum(row_partitions[:i+1]) for i in range(p)]
        col_bounds = [0] + [sum(col_partitions[:i+1]) for i in range(p)]
        
        # Check each block
        for i in range(p):
            for j in range(p):
                # Get block boundaries
                r_start, r_end = row_bounds[i], row_bounds[i+1]
                c_start, c_end = col_bounds[j], col_bounds[j+1]
                
                # Extract block
                block = matrix_L[r_start:r_end, c_start:c_end]
                
                # Check if it's an anti-diagonal block
                is_anti_diagonal = (i + j == p - 1)
                
                # If not anti-diagonal, block must be all zeros
                if not is_anti_diagonal:
                    if np.any(block != 0):
                        return False
                        
        return True
    
    def calculate_objective(self, num_rows_U, row_partitions, p):
        """
        Calculate the objective function value.
        
        Args:
            num_rows_U: number of rows in upper half U
            row_partitions: list of row partition sizes
            p: number of partitions
            
        Returns:
            float: objective value
        """
        # Find anti-diagonal blocks and their row counts
        anti_diagonal_rows = []
        for i in range(p):
            j = p - 1 - i  # anti-diagonal condition: i + j = p - 1
            if j >= 0 and j < p:
                anti_diagonal_rows.append(row_partitions[i])
        
        max_anti_diagonal_rows = max(anti_diagonal_rows) if anti_diagonal_rows else 0
        return num_rows_U + max_anti_diagonal_rows
    
    def try_partition(self, matrix, p_values=[2, 3, 4]):
        """
        Try different values of p and find the best partition.
        
        Args:
            matrix: permuted matrix to partition
            p_values: list of p values to try
            
        Returns:
            tuple: (best_objective, best_solution_info)
        """
        best_obj = float('inf')
        best_info = None
        
        for p in p_values:
            # Try different ways to split rows and columns into p partitions
            for num_rows_U in range(min(self.m, self.m - p + 1)):
                matrix_U = matrix[:num_rows_U, :]
                matrix_L = matrix[num_rows_U:, :]
                
                if matrix_L.shape[0] == 0:
                    continue
                
                # Try different row partitions for L
                rows_L = matrix_L.shape[0]
                cols_L = matrix_L.shape[1]
                
                # Generate possible row partitions
                for row_partitions in self.generate_partitions(rows_L, p):
                    # Generate possible column partitions
                    for col_partitions in self.generate_partitions(cols_L, p):
                        # Check if this partition satisfies the anti-diagonal structure
                        if self.check_anti_diagonal_structure(matrix_L, p, row_partitions, col_partitions):
                            obj = self.calculate_objective(num_rows_U, row_partitions, p)
                            if obj < best_obj:
                                best_obj = obj
                                best_info = {
                                    'p': p,
                                    'num_rows_U': num_rows_U,
                                    'row_partitions': row_partitions,
                                    'col_partitions': col_partitions,
                                    'objective': obj
                                }
        
        return best_obj, best_info
    
    def generate_partitions(self, total, p):
        """
        Generate all ways to partition a number into p positive parts.
        
        Args:
            total: total number to partition
            p: number of parts
            
        Returns:
            list: list of partitions
        """
        if p == 1:
            return [[total]]
        
        partitions = []
        for i in range(1, total - p + 2):
            for rest in self.generate_partitions(total - i, p - 1):
                partitions.append([i] + rest)
        
        return partitions
    
    def check_k_constraint(self, row_perm, num_rows_U, row_partitions):
        """
        Check if the first k rows constraint is satisfied.
        
        Args:
            row_perm: row permutation
            num_rows_U: number of rows in U
            row_partitions: row partition sizes
            
        Returns:
            bool: True if constraint is satisfied
        """
        # Find where the first k rows ended up
        first_k_positions = []
        for i in range(min(self.k, len(row_perm))):
            pos = list(row_perm).index(i)
            first_k_positions.append(pos)
        
        # Check if all first k rows are in U or first row partition of L
        first_partition_end = num_rows_U + (row_partitions[0] if row_partitions else 0)
        
        for pos in first_k_positions:
            if pos >= first_partition_end:
                return False
        
        return True
    
    def solve(self, max_permutations=1000):
        """
        Solve the matrix partitioning problem.
        
        Args:
            max_permutations: maximum number of permutations to try
            
        Returns:
            dict: solution information
        """
        # For large matrices, we'll use a heuristic approach
        # Try a limited number of random permutations
        
        row_indices = list(range(self.m))
        col_indices = list(range(self.n))
        
        tried_permutations = 0
        
        # Try identity permutation first
        permutation_pairs = [(row_indices, col_indices)]
        
        # Add some random permutations
        np.random.seed(42)
        for _ in range(min(max_permutations - 1, 99)):
            row_perm = np.random.permutation(row_indices).tolist()
            col_perm = np.random.permutation(col_indices).tolist()
            permutation_pairs.append((row_perm, col_perm))
        
        for row_perm, col_perm in permutation_pairs:
            tried_permutations += 1
            
            # Apply permutation
            permuted_matrix = self.original_matrix[np.ix_(row_perm, col_perm)]
            
            # Try partitioning with different p values
            obj, info = self.try_partition(permuted_matrix)
            
            if info and obj < self.best_objective:
                # Check k constraint
                if self.check_k_constraint(row_perm, info['num_rows_U'], info['row_partitions']):
                    self.best_objective = obj
                    self.best_solution = {
                        'row_permutation': row_perm,
                        'col_permutation': col_perm,
                        'permuted_matrix': permuted_matrix,
                        **info
                    }
            
            if tried_permutations >= max_permutations:
                break
        
        return self.best_solution
    
    def visualize_solution(self):
        """
        Visualize the best solution found.
        """
        if not self.best_solution:
            print("No solution found.")
            return
        
        sol = self.best_solution
        matrix = sol['permuted_matrix']
        
        print(f"Best solution found:")
        print(f"Objective value: {sol['objective']}")
        print(f"p = {sol['p']}")
        print(f"Number of rows in U: {sol['num_rows_U']}")
        print(f"Row partitions: {sol['row_partitions']}")
        print(f"Column partitions: {sol['col_partitions']}")
        print(f"Row permutation: {sol['row_permutation']}")
        print(f"Column permutation: {sol['col_permutation']}")
        
        print("\nPermuted matrix:")
        print(matrix)
        
        # Show the block structure
        print("\nBlock structure:")
        num_rows_U = sol['num_rows_U']
        if num_rows_U > 0:
            print(f"Upper half U (rows 0-{num_rows_U-1}):")
            print(matrix[:num_rows_U, :])
        
        if num_rows_U < matrix.shape[0]:
            matrix_L = matrix[num_rows_U:, :]
            print(f"Lower half L (rows {num_rows_U}-{matrix.shape[0]-1}):")
            print(matrix_L)
            
            # Show individual blocks
            row_partitions = sol['row_partitions']
            col_partitions = sol['col_partitions']
            p = sol['p']
            
            row_bounds = [0] + [sum(row_partitions[:i+1]) for i in range(p)]
            col_bounds = [0] + [sum(col_partitions[:i+1]) for i in range(p)]
            
            print(f"\nBlocks in lower half L:")
            for i in range(p):
                for j in range(p):
                    r_start, r_end = row_bounds[i], row_bounds[i+1]
                    c_start, c_end = col_bounds[j], col_bounds[j+1]
                    block = matrix_L[r_start:r_end, c_start:c_end]
                    is_anti_diagonal = (i + j == p - 1)
                    print(f"Block ({i},{j}) {'(anti-diagonal)' if is_anti_diagonal else ''}:")
                    print(block)
                    print()

# Example usage
if __name__ == "__main__":
    # Example 1: Small test matrix
    test_matrix = np.array([
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
    
    k = 1  # First k rows must go to U or first row partition of L
    
    solver = MatrixPartitionSolver(test_matrix, k)
    solution = solver.solve(max_permutations=100)
    
    if solution:
        solver.visualize_solution()
    else:
        print("No valid solution found.")