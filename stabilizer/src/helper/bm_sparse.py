import numpy as np
import scipy.sparse as sparse

class BinaryMatrixSparse:
    def __init__(self,bm=None):
        if bm is None:
            self.bm = []
            self.shape = (0,)
        else:
            try:
                self.bm = sparse.csr_matrix(bm)
            except TypeError:
                raise TypeError('The input for BinaryMatrix needs to be a matrix.')
            self.shape = self.bm.shape
        return
    
    def set_matrix(self,bm):
        if isinstance(bm,BinaryMatrixSparse):
            self.bm = bm.bm.copy()
            self.bm = self.bm.tocsr()
            self.shape = self.bm.shape
            return

        try:
            self.bm = sparse.csr_matrix(bm)
        except TypeError:
            raise TypeError('The input for BinaryMatrix needs to be a matrix.')
        self.shape = self.bm.shape
        return
    
    def row_permute(self,indices,trans_mtx_return = False):
        '''
        permute the rows of the matrix based on the indices list given.
        new position: [0, 1, 2, ...]
        old position: [i, j, k, ...]
        old row i -> new row 0

        e.g. if given [1,2,3,0] for a four row dm:
        new row 0 is from old row 3
        new row 1 is from old row 0
        new row 2 is from old row 1
        new row 3 is from old row 2
        
        Input:
            indices:            the indices list for row permutation.
            trans_mtx_return:   If the transformation matrix will be returned.   
        Output:
            transform_mtx:  the transformation matrix M that corresponds to the current row-permutation.
                            the new matrix nbm = M.bm
        '''
        self.bm= self.bm.tocoo()
        row_index = np.array([indices[x] for x in self.bm.row])
        self.bm.row = row_index
        self.bm= self.bm.tocsr()
        if trans_mtx_return:
            elem = [1]*self.shape[0]
            col = range(self.shape[0])
            row = indices
            trans_mtx_return = sparse.coo_matrix((elem,(row,col)))
            return trans_mtx_return
        return
    
    def __add__(self, binary_mtx):
        '''
        addtion in binary with mod 2.
        Input:
            binary_mtx:     the binary matrix for addition.
        Output:
            new_bm:         New BinaryMatrix class instance for the result.
        Exception:
            TypeError:      If the given binary_mtx connot be converted to BinaryMatrix instance.
            IndexError:     If the given binary_mtx is a BinaryMatrix instance, but has different shape than the current matrix.
        '''
        if not isinstance(binary_mtx,BinaryMatrixSparse):
            try:
                binary_mtx = BinaryMatrixSparse(binary_mtx)
            except TypeError:
                raise TypeError('The addition only support types that can be converted to BinaryClass.')
            
        if self.shape != binary_mtx.shape:
            raise IndexError('Addition binary matrices have different dimensions.')
            
        new_bm = self.bm + binary_mtx.bm
        new_bm.data = new_bm.data%2
        return BinaryMatrixSparse(new_bm)

    @staticmethod
    def binary_add(mtx1,mtx2):
        '''
        This function provide a general method for adding two sparse matrix togather in mod 2
        return a binary matrix type instance
        '''
        mtx1 = BinaryMatrixSparse(mtx1)
        mtx3 = (mtx1 + mtx2)
        
        return mtx3

    def __mul__(self,binary_mtx):
        '''
        This provide the matrix dot product with Mod 2 for two binary matrices.
        Note that A * B is A.dot(B) in scipy fasion.
        '''
        if not isinstance(binary_mtx,BinaryMatrixSparse):
            binary_mtx = BinaryMatrixSparse(binary_mtx)
        try:
            temp = self.bm.dot(binary_mtx.bm)
        except ValueError:
            raise ValueError('Dimension mismatch for two binary matrices.')
        temp.data = temp.data%2
        return BinaryMatrixSparse(temp)
    
    def __getitem__(self,i,j=None):
        '''
        Get the row with index i
        Get the element with index (i,j)
        Input:
            i:  row index
            j:  col index (can be ignored to get the whole row)
        Output:
            if only row index is given, it will return a BinaryMatrix instance.
            If both indices are given, it would be in int type.
        '''
        if j is None:
            return BinaryMatrixSparse(self.bm[i])
        else:
            return self.bm[i,j]

    def __str__(self):
        '''
        Print the qubit number.
        Print the Pauli operators for each bits
        '''
        return str(self.bm.toarray())

    def __repr__(self):
        return self.__str__()
    
    def dot(self, new_bm):
        if isinstance(new_bm, BinaryMatrixSparse):
            mtx = self.bm.dot(new_bm.bm)
            mtx.data = mtx.data%2
            return BinaryMatrixSparse(mtx)
        else:
            if sparse.issparse(new_bm):
                mtx = self.bm.dot(new_bm.bm)
                mtx.data = mtx.data%2
                return BinaryMatrixSparse(mtx)

            elif isinstance(new_bm, np.ndarray):
                mtx = self.bm.dot(new_bm)
                mtx = mtx%2
                return BinaryMatrixSparse(mtx)
            else:
                raise TypeError('The type is not supported')

    def row_insert(self,i,j,trans_mtx_return = False):
        '''
        Insert the i-th row to the original j-th row position.
        Need i>j, while i < total rows
        Input:
            i:                  the row index that are waiting for insertion
            j:                  the row index for the insertion place
            trans_mtx_return:   whether the transformation matrix is returned.
        Output:
            The function changes self.tableau
            trans_mtx:          The transformation matrix if trans_mtx_return is True
        Exception:
            ValueError: if the row index is invalid
        '''
        total_row = self.bm.shape[0]
        
        if i>=total_row:
            raise ValueError("The inserted row index is out of range.")
        if i<j:
            raise ValueError("The inserted row index is before the insertion place.")
        
        if i==j:
            if trans_mtx_return:
                return sparse.eye(total_row)
            return
        
        ind1 = list(range(0,j)) + list(range(j+1,i+1))
        ind1.append(j)
        ind1 = ind1 + list(range(i+1,total_row))
        #print(ind1)
        trans_mtx = self.row_permute(ind1,trans_mtx_return)
        
        if trans_mtx_return:
            return trans_mtx
        else:
            return

    def search_nonzero_rows(self,i,start_row = 0, final_row = -1):
        '''
        This function is to select the rows that are nonzero for a specific column index i
        return the row indices that are zero and row indices that are nonzero
        The slice range [start_row:final_row] (the final row index is excluded.)
        Input:
            i:          the column index for seaching nonzero rows
            start_row:  the starting rows to search, default 0
            final_row:  the final rows that finish the search, default -1
        Output:
            [zeros, ones]
            zeros:      the indices for the rows that column i element is 0
            ones:       the indices for the rows that column i element is 1
        '''
        total_col = self.bm.shape[1]
        
        if final_row ==-1:
            final_row = self.bm.shape[0]
        
        if i >= total_col:
            raise ValueError("The column index is out of range.")
        
        column = np.array(self.bm.tocsc().getcol(i).todense()).flatten()
        
        column = column[start_row:final_row]
        
        #print(column)
        
        zeros = np.where(column == 0)[0] + start_row
        ones = np.where(column == 1)[0] + start_row
        return [zeros, ones]
        
    def _reduce_down(self,column_indices,trans_mtx_return=False):
        '''
        This is a helper function for the reduction function.
        The given input column_indices are the columns we are going to calcualtion.
        This helper function use the first row, to add to the other rows.
        '''
        self_copy = BinaryMatrixSparse(self.bm)
        self_copy.bm = self_copy.bm.tolil()
        
        j0 = column_indices[0]
        
        for j in column_indices[1:]:
            ########
            # DEBUG CODE
            #print('j,j0',j,j0)
            #print('mtx[j] = ',mtx[j].todense())
            #print('mtx[j0] = ',mtx[j0].todense())
            ########
            ## res = Tableau.tab_add(mtx[j0],mtx[j])
            res = self_copy[j0] + self_copy[j]
            if res.bm.count_nonzero() == 0:
                self_copy.bm[j,:]=0
            else:
                self_copy.bm[j] = res.bm.copy()
            ########
            # DEBUG CODE
            #print('tab_add = ',Tableau.tab_add(mtx[j0],mtx[j]),' type:',type(Tableau.tab_add(mtx[j0],mtx[j])))
            #print('mtx[j+j0] = ',mtx[j].todense())
        
        #print('mtx = ',mtx.todense())
        #########
        self.set_matrix(self_copy)
        #print('self = ',self.tableau.todense())
        self.bm.eliminate_zeros()

        if trans_mtx_return:
            elem = [1] * (self.shape[0] + len(column_indices) - 1)
            row = list(range(self.shape[0])) + list(column_indices[1:])
            col = list(range(self.shape[0])) + [column_indices[0]] * (len(column_indices) -1)
            trans_mtx = sparse.csr_matrix((elem,(row,col)))
            return trans_mtx

        return
    
    def _reduce_up(self,column_indices,trans_mtx_return=False):
        '''
        This is a helper function for the reduction function.
        The given input column_indices are the columns we are going to calcualtion.
        This helper function use the last row, to add to the other rows.
        '''
        self_copy = BinaryMatrixSparse(self.bm)
        self_copy.bm = self_copy.bm.tolil()

        j0 = column_indices[-1]
        
        for j in column_indices[:-1]:
            res = self_copy[j0] + self_copy[j]
            if res.bm.count_nonzero() == 0:
                self_copy.bm[j,:]=0
            else:
                self_copy.bm[j] = res.bm.copy()
        
        self.set_matrix(self_copy.bm)
        self.bm.eliminate_zeros()
        
        if trans_mtx_return:
            elem = [1] * (self.shape[0] + len(column_indices) - 1)
            row = list(range(self.shape[0])) + list(column_indices[:-1])
            col = list(range(self.shape[0])) + [column_indices[-1]] * (len(column_indices) -1)
            trans_mtx = sparse.csr_matrix((elem,(row,col)))
            return trans_mtx

        return
        

    def reduction(self,trans_mtx_return = False):
        '''
        This function is to perform the Gaussian reduction to the matrix to RREF.
        Here we do not allow column exchange.
        
        Input:
            None
        Output:
            None: the function will change the self.bm
        '''
        row_pointer = 0      ## current row position
        column_pointer = 0   ## current column position
        trans_mtx_list = []
        
        for column_pointer in range(0,self.shape[1]):
            #print(row_pointer)
            _,ones = self.search_nonzero_rows(column_pointer,start_row = row_pointer)
            ## group the rest of the rows by either the column == 1 or not
            if len(ones) == 0:
                ## there is no rows that has "1" on "column_pointer" column position
                ## we skip this column and do nothing
                continue
            else: 
                #print('reduce down ones:',ones)
                ## make sure the row: ones[0] is the only row that have a 1 on the "column_pointer" column position
                trans_mtx1 = self._reduce_down(ones, trans_mtx_return)
                ## insert this row to the current row_pointer pointed rows
                trans_mtx2 = self.row_insert(ones[0],row_pointer,trans_mtx_return)                   
                
                if trans_mtx_return:
                    trans_mtx_list.append(trans_mtx1) 
                    trans_mtx_list.append(trans_mtx2)

                row_pointer += 1
                
                #print(self.tableau.todense())
                
                ## following: we try to use this new row, to elimitate the previous rows that has one on this working column
                _, ones_up = self.search_nonzero_rows(column_pointer,final_row = row_pointer)
                trans_mtx1 = self._reduce_up(ones_up,trans_mtx_return)

                if trans_mtx_return:
                    trans_mtx_list.append(trans_mtx1)
                
                #print(self.tableau.todense())
                
                self.bm.eliminate_zeros()
                
                ## if we are working on the last row, we are done
                if row_pointer == self.shape[0]:
                    break
        
        # self._qubit_relabel(doc_file)
        if trans_mtx_return:
            return trans_mtx_list
        return

    def col_permute(self,indices,trans_mtx_return = False):
        '''
        permute the qubits (relabel the qubits) of the tableau based on the indices list given.
        new position: [0, 1, 2, ...]
        old position: [i, j, k, ...]
        old qubit i -> new qubit 0

        Input:
            indices: the new qubit index for qubit permuting
        Output:
            None, will change the current instance
        
        Note:
            (1) the index list will be checked. If the list length does not match the qubit number, a TypeError will be raised.
        '''
        if len(indices) != self.shape[1]:
            raise TypeError('The permutation list does not match the column number.')

        self.bm = self.bm.tocoo()
        col_index = np.array([indices[x] for x in self.bm.col])
        self.bm.col = col_index
        self.bm = self.bm.tocsr()
        if trans_mtx_return:
            elem = [1]*self.shape[1]
            row = range(self.shape[1])
            col = indices
            trans_mtx = sparse.coo_matrix((elem,(row,col)))
            return trans_mtx
        return
    
    def _col_relabel(self):
        '''
        This function is to allow permutation on the columns to make the pivot elements all on the leftmost positions
        [[I  A ]
         [0  0 ]]
        where I has the maximum dimension.
        '''
        col_index = np.zeros(self.shape[1])
        bm_mtx = self.bm.copy()

        filled_set = set({})
        for j in range(0,self.shape[0]):
            col = bm_mtx[j].nonzero()[1]
            # print('j',j)
            if len(col) == 0:
                j -= 1
                break
            else:
                col_index[col[0]] = j
                filled_set.add(col[0])
        
        # print('index: ',qubit_index)

        for i in range(0,self.shape[1]):
            if i in filled_set:
                continue
            else:
                # print('j',j)
                j+=1
                col_index[i] = j
        
        # print('index: ',col_index)
        col_index = col_index.astype(int)
        trans_mtx = self.col_permute(col_index,True)

        return trans_mtx
    
    def non_zero_lines(self,reduced = True):
        '''
        this function solves the number of nonzero lines of the binary matrix.
        '''
        if not reduced: 
            self.reduction()
        
        self.bm.eliminate_zeros()

        row_index = self.shape[0]-1

        while row_index >= 0 and len(self.bm[row_index].data) == 0:
            row_index -= 1        
        return row_index + 1

    def rank(self,reduced = True):
        return self.non_zero_lines(reduced)
    
    def nullspace(self):
        '''
        This function solves the nullspace of the binary matrix.
        The matrix will be at first reduce to rref form. The corresponding transformation matrix is recorded.
        The nullspace of the transformed matrix can be solved.
        Then the nullspace will be converted back.

        Input:
            None
        Output:
            bm_null:    This is a BinaryMatrix instance. Each row is a vector in the null-space.
        '''
        self.reduction()
        ## the number of nonzero lines in the matrix
        r = self.non_zero_lines()

        ## if the nonzero lines are equal to the second shape of the binary matrix
        ## there is no non-trivial nullspace.
        ## The result null-space will be a single zero vector

        if r == self.shape[1]:
            return BinaryMatrixSparse(np.zeros([self.shape[1],1]).astype(int))

        col_trans_mtx = self._col_relabel()
        bm_mtx = self.bm.tolil()

        a_mtx = bm_mtx[:r,r:]

        ## for debug
        #print('a_mtx shape: '+str(a_mtx.shape))

        i_mtx = sparse.eye(self.shape[1] - r)

        null_space_mtx = sparse.bmat([[a_mtx],[i_mtx]])
        null_space_mtx = col_trans_mtx.dot(null_space_mtx)

        bm_null = BinaryMatrixSparse(null_space_mtx.T)
        self.bm = self.bm.dot(col_trans_mtx.T) ## convert the bm back before relabel the columns
        return bm_null, (col_trans_mtx)

    def transpose(self):
        trans = self.bm.T
        trans = trans.astype(int)
        temp = BinaryMatrixSparse(trans) 
        return temp

    def inverse(self):
        '''
        find the inverse of the current binary matrix
        '''

        if self.shape[0] != self.shape[1]:
            raise ValueError('This matrix is not square.')

        length = self.shape[0]

        bm_expanded = np.zeros((length,2*length), dtype=int)

        bm_expanded[:, :length] += self.bm
        bm_expanded[:, length:] += np.eye(length, dtype=int)
        
        bm_new_for_expand = BinaryMatrixSparse(bm_expanded)
        bm_new_for_expand.reduction()

        ## check if the original matrix becomes identity matrix or not
        if not np.array_equal(bm_new_for_expand.bm[:, :length].toarray(), np.eye(length, dtype=int)):
            print(bm_new_for_expand[:, :length])
            raise ValueError('The matrix is not invertable')
        
        else:
            return bm_new_for_expand[:, length:]