import numpy as np
import scipy.sparse as sparse

#from Stabilizer.BinaryMatrixSparse import BinaryMatrixSparse as BM
from .bm_sparse import BinaryMatrixSparse as BM

#OPERATOR_CONST = ['I','Z','X','Y']
OPERATOR_CONST = {0:'I',1:'Z',2:'X',3:'Y','I':0,'Z':1,'X':2,'Y':3,'i':0,'z':1,'x':2,'y':3}
ANTI_SYM_TENSOR = {(1,2):0, (2,3):0,(3,1):0,(2,1):1,(3,2):1,(1,3):1}
SINGLE_GATE_SET = set({'X','x','Y','y','Z','z','H','h','S','s'})
TWO_GATE_SET = set({'cnot','CNOT','CZ','cz','SWAP','swap'})

class Tableau:
    def __init__(self,tab):
        ## dimensional check:
        ##   1). the last bit should be to track the phase of the corresponding stabilizer operator, so column number must be odd
        ##   2). the row number should <= column number // 2
        ## 
        self.tableau = sparse.csr_matrix(tab)
        row_number, col_number = self.tableau.shape
        
        if col_number%2 == 0:
            raise ValueError('The input tableau dimension is incorrect. Each row should be odd.')
        if row_number > col_number //2:
            raise ValueError('The input tableau has more rows than the qubit number.')
        
        self.qubits = col_number // 2
        self.stab_counts = row_number
        
        self.shape = self.tableau.shape
        
        return
    
    def tab(self):
        '''
        This function gives the tableau matrix with the phase bits
        Return: 
            sparse matrix for the tableau (with the phase bit)
            Dimension [s x (2n + 1)], where s -> stabilizer counts, n -> qubit number
        '''
        return self.tableau
    
    def stabilizers(self):
        '''
        This function returns the tableau for the stabilizers, without the phase
        Return:
            Sparse matrix for the tableau (without phase bit)
            The sparse matrix is in coo_matrix format.
            Dimension [s x (2 n)], where s -> stabilizer counts, n -> qubit number
        '''
        return self.tableau[:,:-1]
    
    def phase(self):
        '''
        This function returns the phase bit for each stabilizers
        Return:
            Sparse matrix in coo_matrix format
            Dimension: [s x 1], where s -> stabilizer counts
        '''
        return self.tableau[:,-1]
    
    def get_index(self,row = None, input_stab = None):
        '''
        This function calculate the index matrix for the stabilizer,
        we assign each single-qubit Pauli operator an index, which is given in the order of OPERATOR_CONST
        Input:
            row: the row number in the current instance that we want to convert to index representation
            input_stab: the input stabilizer that are going to be converted
            
            if no input is given, the function calcualte the index representation of the tableau of the current instance
            if both row and input_stab are given, only row input will be processed.
        '''
        if row is None and input_stab is None:
            stab = self.stabilizers()
            x_part = stab[:,:self.qubits]
            z_part = stab[:,self.qubits:]
            
            index_mtx = x_part * 2 + z_part
            return index_mtx
        elif input_stab is None:
            try:
                row = int(row)
            except:
                raise TypeError("get_index row number needs to be integer.")
            
            stab = self.tableau[row,:self.qubits] * 2 + self.tableau[row,self.qubits:-1]
            return stab
        else:
            in_stab = sparse.csr_matrix(input_stab)
            if in_stab.shape[1] % 2 == 1:
                print("Input stab with phase bit")
                in_stab = in_stab[:,:-1]
            else:
                print("Input stab without phase bit")
            
            bits = in_stab.shape[1] // 2
            
            x_part = in_stab[:,:bits]
            z_part = in_stab[:,bits:]
            
            index_mtx = x_part * 2 + z_part
            return index_mtx
        
    def readout(self,row):
        '''
        This function gives a readout output (+/- Pauli operators for each bits) of a certain row of the stabilier 
        Input:
            row: the integer number for the row that want to read out in the tableau
        Output:
            str_out: a string that contains the sign and Pauli operators for each bits
        '''
        try:
            row = int(row)
        except:
            raise TypeError('readout need an integer row number.')
        row_tab_index = self.tableau[row,:self.qubits]*2 + self.tableau[row,self.qubits:-1]
        row_tab_index = row_tab_index.toarray().flatten()
        
        row_tab_index_str = [OPERATOR_CONST[x] for x in row_tab_index]
        
        row_sign = self.tableau[row,-1]
        if row_sign == 0:
            str_out = '+'
        else:
            str_out = '-'
        
        str_out += ''.join(row_tab_index_str)
        return str_out
    
    def __str__(self):
        '''
        Print the qubit number.
        Print the Pauli operators for each bits
        '''
        row = self.tableau.shape[0]
        str_out = '-'*max(self.qubits,30)+'\n'
        str_out += 'Qubit number: '+str(self.qubits)+'\n'
        str_out += 'Stabilizer count: '+str(row) + '\n'
        str_out += 'Stabilizers: \n'
        
        for i in range(0,row):
            str_out += self.readout(i)
            str_out += '\n'
        str_out += '-'*max(self.qubits,30)+'\n'
        
        return str_out

    def __repr__(self):
        return self.__str__()
    
    def __getitem__(self,i,j=None):
        if j is None:
            return self.tableau[i]
        else:
            return self.tableau[i,j]
    
    @staticmethod
    def convert_back(string_in):
        '''
        This function convert a string of Pauli Gates, to stabilizer formalism
        Input:
            string_in:  The string of Pauli Gates to convert to the tableau formalism
        Output:
            out_stab:   The stabilizer for output
        Note:
            (1) The input string can be in upper or lower letters
            (2) The string can be start from + or -, which will give the phase bit
            (3) If the string does not start by the sign, we will assume the sign is +
            (4) The input can also be a list of string. The output will be a list of Tableau instances for each line
            (5) If the input is neither a string, nor a list of string (list or ndarray), a TypeError will be raised.
        '''
        if isinstance(string_in,str):
            ## consider the first char, if it is a sign character
            phase_bit = 0
            gates = string_in
            if string_in[0] == '+' or string_in[0] == '-':
                gates = string_in[1:]
                if string_in[0] == '-':
                    phase_bit = 1
            try:
                index_list = [OPERATOR_CONST[x] for x in gates]
            except KeyError:
                raise KeyError('Input string needs to be only { X, Y, Z, I }.')
            
            index_list = np.array(index_list)
            tab = np.zeros(len(index_list)*2+1)
            tab[-1] = phase_bit
            tab[:len(index_list)] = index_list//2
            tab[len(index_list):-1] = index_list%2
            tab = tab.astype(int)

            out_tab = Tableau(tab)
            return out_tab
        elif isinstance(string_in,list) or isinstance(string_in,np.ndarray):
            ## Here we assume the input 
            out_tab_list = []
            for string in string_in:
                out_tab_list.append(Tableau.convert_back(string))
            
            return out_tab_list
        else:
            raise TypeError('The input type can only be string, list, ndarray.')

    @staticmethod
    def binary_add(mtx1,mtx2):
        '''
        This function provide a general method for adding two sparse matrix togather in mod 2
        '''
        mtx1 = sparse.csr_matrix(mtx1)
        mtx2 = sparse.csr_matrix(mtx2)
        
        mtx3 = mtx1 + mtx2
        mtx3.data = mtx3.data % 2
        mtx3.eliminate_zeros()
        
        return mtx3
    
    @staticmethod
    def tab_add(line1,line2,index_1 = None, index_2 = None, output_format='s'):
        '''
        This function provide a general method for adding two stabilizers' tableau representation (with phase bit) togather
        Inputs:
            line1, line2: two lines of stabilizers (with phase)
            index_1, index_2 : The index vector for the given two lines
                               These two inputs can be omitted.
            output_format: 's':  for sparse array stabilizer expression
                           'n': for normal (dense list)
            
            Warning: 1). If the index_1 and index_2 are given, the correctness of these two inputs are not checked by the function!
                         So be sure you are using the authentic index matrices.
                     2). There is no dimension check in the function.
                         We automatically assume that both the inputs are in the stabilizer (with phase) convension.
                         If the given dimension does not match (2nd dimension is not odd), it may cause unexpected results
                     3). We assume both inputs are in sparse matrix type
        '''
        line1 = sparse.csr_matrix(line1)
        line2 = sparse.csr_matrix(line2)
        
        bits = line1.shape[1] // 2
        
        stab1 = line1[0,:-1]
        stab2 = line2[0,:-1]
        
        if index_1 is None:
            index_1 = line1[0,:bits] * 2 + line1[0,bits:-1]
        if index_2 is None:
            index_2 = line2[0,:bits] * 2 + line2[0,bits:-1]

        phase1 = line1[0,-1]
        phase2 = line2[0,-1]
        
        stab3 = stab1 + stab2
        stab3.data = stab3.data % 2
        stab3.eliminate_zeros()
        
        index_3 = (stab3[0,:bits]*2 + stab3[0,bits:])
        
        col1 = index_1.nonzero()[1]  ## the column number for stab1 (self) that are not I
        col2 = index_2.nonzero()[1]  ## the column number for stab_in that are not I
        col3 = index_3.nonzero()[1]  ## the column number for result stabilizer that are not I (self and input are identical Pauli ops)
        
        intercect = np.intersect1d(np.intersect1d(col1,col2),col3)  ## these are the only bits that contributes to overall phase
        ## each bit of the previous intercect list at least contributes to an i 
        if len(intercect) % 2 != 0:
            print("Warning: the two stabilizers do not commute.")
        
        ## Here we want to determine the extra minus sign for the multiplication process
        ## 0 for + and 1 for -
        sign_list = [ANTI_SYM_TENSOR[(index_1[0,x],index_2[0,x])] for x in intercect]
        
        #print(phase1,phase2,col3,len(col3)//2,sign_list,np.sum(sign_list))
        
        new_sign = (phase1 + phase2 + len(intercect)//2 + np.sum(sign_list))%2
        
        #print("new signl;", new_sign)
        
    
        #new_stab_vect = sparse.bmat([[stab3,sparse.csr_matrix(new_sign)]])
        # print('debug: stab3')
        # print(new_sign)
        # print(type(new_sign))
        if not sparse.issparse(sign_list):
            # print('here')
            new_stab_vect = sparse.bmat([[stab3,[new_sign]]])
        else:
            # print('there')
            new_stab_vect = sparse.bmat([[stab3,(new_sign)]])
        
        if output_format == 's':
            return new_stab_vect.tocsr()
        elif output_format == 'n':
            return new_stab_vect.todense()
    
    @staticmethod
    def commute(line1,line2):
        '''
        This function tests the commutation relation of two stabilizers given.
        Here we assume that the input are the tableau representation (with phase)
        In the function, we test the size of the two lines, if they did not match, a ValueError will be raised.
        Input:
            line1, line2: input tableau of two stabilizers with phase bits
        Output:
            Boolean type
        '''
        
        line1 = sparse.csr_matrix(line1)
        line2 = sparse.csr_matrix(line2)
        
        ## test length of the two tableau
        if line1.shape[1] != line2.shape[1]:
            raise ValueError("The input stabilizers have different qubit size.")
        
        stab1 = line1[0,:-1]
        stab2 = line2[0,:-1]
        
        bits = stab1.shape[1] // 2
        
        imtx = sparse.eye(bits)
        gamma_mtx = sparse.bmat([[None,imtx],[imtx,None]])
        
        res = stab1.dot(gamma_mtx).dot(stab2.T)[0,0]
        res = res % 2
        if res == 0:
            return True
        else:
            return False
        
    def row_permute(self,indices):
        '''
        permute the rows of the tableau based on the indices list given.
        new position: [0, 1, 2, ...]
        old position: [i, j, k, ...]
        old row i -> new row 0
        
        Input:
            indices:    the indices list for row permutation.
        Output:
            None:       the function changes the self.tableau
        '''
        self.tableau = self.tableau.tocoo()
        row_index = np.array([indices[x] for x in self.tableau.row])
        #print(self.tableau.row)
        #print(row_index)
        self.tableau.row = row_index
        #print(self.tableau.row)
        self.tableau = self.tableau.tocsr()
        return
    
    def row_insert(self,i,j,doc_file = None):
        '''
        Insert the i-th row to the original j-th row position.
        Need i>j, while i < total rows
        Input:
            i:          the row index that are waiting for insertion
            j:          the row index for the insertion place
            doc_file:   the recording document file name.
        Output:
            None:       The function changes self.tableau
        Exception:
            ValueError: if the row index is invalid
        '''
        total_row = self.tableau.shape[0]
        
        if i>=total_row:
            raise ValueError("The inserted row index is out of range.")
        if i<j:
            raise ValueError("The inserted row index is before the insertion place.")
        
        if i==j:
            return
        
        ind1 = list(range(0,j)) + list(range(j+1,i+1))
        ind1.append(j)
        ind1 = ind1 + list(range(i+1,total_row))
        #print(ind1)
        self.row_permute(ind1)
        
        if doc_file is None:
            return
        else:
            with open(doc_file,'a') as fout:
                fout.write('Line insertion:\t'+str(i)+'\t'+str(j)+'\n')
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
        total_col = self.tableau.shape[1]
        
        if final_row ==-1:
            final_row = self.tableau.shape[0]
        
        if i >= total_col:
            raise ValueError("The column index is out of range.")
        
        column = np.array(self.tableau.tocsc().getcol(i).todense()).flatten()
        
        column = column[start_row:final_row]
        
        #print(column)
        
        zeros = np.where(column == 0)[0] + start_row
        ones = np.where(column == 1)[0] + start_row
        return [zeros, ones]
        
    def _reduce_down(self,column_indices,doc_file = None):
        '''
        This is a helper function for the reduction function.
        The given input column_indices are the columns we are going to calcualtion.
        This helper function use the first row, to add to the other rows.
        '''
        mtx = self.tableau.tolil()
        
        j0 = column_indices[0]
        
        for j in column_indices[1:]:
            ########
            # DEBUG CODE
            #print('j,j0',j,j0)
            #print('mtx[j] = ',mtx[j].todense())
            #print('mtx[j0] = ',mtx[j0].todense())
            ########
            res = Tableau.tab_add(mtx[j0],mtx[j])
            if res.count_nonzero() == 0:
                mtx[j,:]=0
            else:
                mtx[j] = res.copy()
            ########
            # DEBUG CODE
            #print('tab_add = ',Tableau.tab_add(mtx[j0],mtx[j]),' type:',type(Tableau.tab_add(mtx[j0],mtx[j])))
            #print('mtx[j+j0] = ',mtx[j].todense())
        
        #print('mtx = ',mtx.todense())
        #########
        self.tableau = mtx.tocsr()
        #print('self = ',self.tableau.todense())
        self.tableau.eliminate_zeros()
        if doc_file is None:
            return
        else:
            with open(doc_file,'a') as fout:
                col_write = ' '.join([str(x) for x in column_indices[1:]])
                fout.write('Reduce down:\t'+str(j0)+'\t'+col_write+'\n')
            return
        
    
    def _reduce_up(self,column_indices,doc_file = None):
        '''
        This is a helper function for the reduction function.
        The given input column_indices are the columns we are going to calcualtion.
        This helper function use the last row, to add to the other rows.
        '''
        mtx = self.tableau.tolil()
        
        j0 = column_indices[-1]
        
        for j in column_indices[:-1]:
            res = Tableau.tab_add(mtx[j0],mtx[j])
            if res.count_nonzero() == 0:
                mtx[j,:]=0
            else:
                mtx[j] = res.copy()
        
        self.tableau = mtx.tocsr()
        self.tableau.eliminate_zeros()
        
        if doc_file is None:
            return
        else:
            with open(doc_file,'a') as fout:
                col_write = ' '.join([str(x) for x in column_indices[:-1]])
                fout.write('Reduce up:\t'+str(j0)+'\t'+col_write+'\n')
            return
        
    def reduction(self,doc_file = None, qubit_exchange = True):
        '''
        This function is to perform the Gaussian reduction to the matrix.
        Here we do allow column exchange.
        
        Do the Gaussian reduction on the tableau matrix
        We will try to simplify the matrix in the form:
            [[I, A | B, C],
             [0, 0 | D, E]]
            where I is an Identity matrix.
        
        Input:
            None
        Output:
            None: the function will change the self.tableau
        '''
        row_pointer = 0      ## current row position
        column_pointer = 0   ## current column position
        
        for column_pointer in range(0,self.tableau.shape[1]):
            #print(row_pointer)
            _,ones = self.search_nonzero_rows(column_pointer,start_row = row_pointer)
            ## group the rest of the rows by either the column == 1 or not
            if len(ones) == 0:
                ## there is no rows that has "1" on "column_pointer" column position
                ## we skip this column and do nothing
                continue
            else: 
                #print('reduce down ones:',ones)
                self._reduce_down(ones, doc_file)    ## make sure the row: ones[0] is the only row that have a 1 on the "column_pointer" column position
                self.row_insert(ones[0],row_pointer,doc_file)  ## insert this row to the current row_pointer pointed rows

                row_pointer += 1
                
                #print(self.tableau.todense())
                
                ## following: we try to use this new row, to elimitate the previous rows that has one on this working column
                _, ones_up = self.search_nonzero_rows(column_pointer,final_row = row_pointer)
                self._reduce_up(ones_up,doc_file)
                
                #print(self.tableau.todense())
                
                self.tableau.eliminate_zeros()
                
                ## if we are working on the last row, we are done
                if row_pointer == self.tableau.shape[0]:
                    break
        
        if qubit_exchange:
            self._qubit_relabel(doc_file)
        return

    def merge(self,new_stab):
        '''
        Merge two systems into one. The new system is initilized by the stabilizer (1D) "new_stab"
        The dimension of the input stabilizer will be checked.
            If the input stabilzer has more than one line (stabilizer_count > 1), a ValueError will be raised. 
        We require the input new_stab is an instance of Tableau group.
            If not, it will be converted to a Tableau class.
        
        We add the input tableau to the end of the system
        e.g., the current tableau [x|z], the new stabilizer [x1|z1]
            the resulting tableau after the merge [x x1 | z z1 | sign]
        Input:
            new_stab: the new stabilizers which determine the state for the new qubits to be added in
        Output:
            new_matrix: a new Tableau instance include the merged new qubits.
        Exception:
            ValueError: If the input Tableau has more than one stabilizers.
        '''
        stabs = self.stabilizers()
        stabs = stabs.tolil()
        stabs_x = stabs[:,:self.qubits]
        stabs_z = stabs[:,self.qubits:]
        
        phase_bits = self.phase()
        
        ## Check if the input stabilizer is in Tableau class, if not, convert it.
        if not isinstance(new_stab,Tableau):
            new_stab = Tableau(new_stab)
        
        ## Check the number of the input stabilizers,
        ## if more than one, raise a ValueError
        if new_stab.stab_counts > 1:
            raise ValueError("The merge function cannot more than one stabilizers as the input.")

        new_stab_vec = new_stab.stabilizers()
        new_stab_vec = new_stab_vec.tolil()
        
        new_stab_x = new_stab_vec[0,:new_stab.qubits]
        new_stab_z = new_stab_vec[0,new_stab.qubits:]
        
        const_i = np.full((self.stab_counts),1)
        
        #print(new_stab_vec,const_i,new_stab_x.todense())
        
        new_x = sparse.kron(const_i,new_stab_x).reshape(self.stab_counts,new_stab.qubits)
        #new_x.data = new_x.data%2
        
        new_z = sparse.kron(const_i,new_stab_z).reshape(self.stab_counts,new_stab.qubits)
        #new_z.data = new_z.data%2
        
        new_phase_bits = sparse.coo_matrix((phase_bits.todense() + new_stab.phase()[0,0]))
        new_phase_bits.data = new_phase_bits.data%2

        #sparse.bmat([[stabs_x,new_x,stabs_z,new_z]])
        
        new_matrix = sparse.bmat([[stabs_x,new_x,stabs_z,new_z,new_phase_bits]])
        new_matrix.eliminate_zeros()
        new_matrix = new_matrix.astype(int)
        
        return Tableau(new_matrix)
    
    def is_commute(self,stab_in):
        '''
        This function calculate whether the given stabilizer(s) commute with the current tableau 
        Will give a list of binary values,  0 -> commute
                                            1 -> anti-commute
        The first dimension corresponds to the input stabilizer(s)
        The second dimension corresponds to the self stabilizers
                                            
        Input:
            stab_in: this is the stabilizer(s) input, which is to determine the commutation relation with the current tableau
                If the stab_in is a Tableau class instance, we directly calculate the commutation relation
                If the stab_in is not a Tableau instance, it will be converted into Tableau instance at first.
        Output:
            commute_relation: the resulting commutation relation, 0 -> commute, 1 -> anti-commute
        Note:
            If the two tableaus have different qubit numbers, the smaller tableau will be extended by appending extra I's.
            A warning will be printed to terminal.
        '''
        
        if not isinstance(stab_in,Tableau):
            stab_in = Tableau(stab_in)
        
        if stab_in.qubits != self.qubits:
            print("Warning: The two stablizers do not have the same qubit numbers.")
            qubit_diff = abs(self.qubits - stab_in.qubits)
            const_0 = np.full(qubit_diff*2 + 1, 0)          ## extra I's
            if stab_in.qubits < self.qubits:
                ## The input stab has smaller qubit number
                ## We need to fill in extra I's to the input stabilizer
                stab_in_new = stab_in.merge(const_0)            ## fill in the extra I's to stab_in
                stab2 = stab_in_new.stabilizers()               ## get the stablizers
                stab1 = self.stabilizers()
                qubit_num = self.qubits
            else:
                ## The input stabilizers have larger qubit number
                ## need to fill in the self stabilizer
                self_new = self.merge(const_0)                  ## fill in the extra I's to stab_in
                stab1 = self_new.stabilizers()                  ## get the stablizers for the new self instance
                stab2 = stab_in.stabilizers()
                qubit_num = stab_in.qubits
        else:
            ## the two stabilizers has the same number of qubits
            stab1 = self.stabilizers()
            stab2 = stab_in.stabilizers()
            qubit_num = self.qubits
        
        ## Construct the intertwine matrix (gamma) for inner product
        ## gamma = [[0 | I],
        ##          [I | 0]]
        qubit_I_mtx = sparse.eye(qubit_num)
        gamma = sparse.bmat([[None,qubit_I_mtx],[qubit_I_mtx,None]]).tocsr()

        commute_relation = sparse.coo_matrix(stab2.dot(gamma).dot(stab1.T)).astype(int)
        commute_relation.data = commute_relation.data % 2 
        return commute_relation

    def append(self, stab_in, doc_file = None):
        '''
        This function append a new stabilizer into the tableau.
        The new stabilizer(s) needs to commute with all the stabilizers inside the system
        The new stabilizers are required to have the same qubit number as the tableau qubit number
        Input:
            stab_in:    The stabilizer(s) that will be appended into the Tableau
        Output:
            None:       self.tableau will be changed
        Exception:
            ValueError: if the input stabilizer does not have the same qubit number
                        if the input stab_in cannot be converted into Tableau class instance
        '''
        ## check the new line dimensions:
        if not isinstance(stab_in,Tableau):
            stab_in = Tableau(stab_in)
        
        ## Check the qubit number of the input stabilizer(s)
        if stab_in.qubits != self.qubits:
            raise ValueError("Tableau.append function requires the input stabilizers have the same qubit numbers.")
        
        ## deal the case that there is only one stabilizer as the input
        commutation_stabilizers = []
        comm_pos = []
        for i in range(0,stab_in.stab_counts):
            commutation = self.is_commute(stab_in[i])      ## calculate the commutation relation between the current insert line (measurement) 
                                                        ## with all the stabilizer in the Tableau
            commutation = commutation.toarray().flatten()

            if not any(commutation.astype(bool)):
                ## this means it commutes with all the stabilizers inside the Tableau
                ## we will append the stabilizer into the commutation_stabilizers list
                commutation_stabilizers.append(stab_in[i].tocoo())
                comm_pos.append(i)
        
        ## append all the commutated stablizers (inside commutation_stabilizers list)
        mtx = self.tableau.tocoo()
        row = list(mtx.row)
        col = list(mtx.col)
        elem = list(mtx.data)

        row0 = self.stab_counts
        for i in range(0,len(commutation_stabilizers)):
            row += [i+row0]*commutation_stabilizers[i].count_nonzero()
            col += list(commutation_stabilizers[i].col)
            elem += list(commutation_stabilizers[i].data)
        
        size = (len(commutation_stabilizers)+row0,self.qubits*2+1)

        new_mtx = sparse.csr_matrix((elem,(row,col)), shape=size)
        self.tableau = new_mtx
        self.stab_counts = size[0]

        if not (doc_file is None):                       ## record the action
            with open(doc_file,'a') as fout:
                output_text = 'Append Action: \n'
                for i in range(0,len(comm_pos)):
                    output_text += 'Append: ' + '\t' + stab_in.readout(i) + '\n'
                fout.write(output_text)
        return

    def measure(self,stab_in, doc_file = None):
        '''
        This function try to measure a new stabilizer (or several stablizers) to the tableau states.
        Input:
            stab_in: the stabilizer(s) that will be measured.
        Output:
            None: the self.tableau will be changed.
        Note:
            (1) The input stabilizers needs to be in Tableau class.
                If not, it will be converted to Tableau class at first
            (2) The input stabilizers should have the same qubit number as the current Tableau.
                If the qubit number is different, a ValueError will be raised.
        '''
        ## check the new line dimensions:
        if not isinstance(stab_in,Tableau):
            stab_in = Tableau(stab_in)
        
        ## Check the qubit number of the input stabilizer(s)
        if stab_in.qubits != self.qubits:
            raise ValueError("Tableau.append function requires the input stabilizers have the same qubit numbers.")

        ## deal the case that there is only one stabilizer as the input
        if stab_in.stab_counts == 1:
            commutation = self.is_commute(stab_in)      ## calculate the commutation relation between the current insert line (measurement) 
                                                        ## with all the stabilizer in the Tableau
            commutation = commutation.toarray().flatten()

            if not any(commutation.astype(bool)):
                ## this means it commutes with all the stabilizers inside the Tableau
                return
            else:
                ## this means there are stabilizers that the not commute with the input stabilizer
                if np.sum(commutation) == 1:
                    position = np.where(commutation ==1)[0][0]       ## The position for the anti-commute stabilizer in the Tableau
                    self.tableau = self.tableau.tolil()
                    self.tableau[position] = stab_in[0].copy()      ## copy the input stabilizer into the place that are not commute with the input one 
                    
                    if not(doc_file is None):                       ## record the action
                        with open(doc_file,'a') as fout:
                            output_text = 'Measure Position: '+'\t' + str(position) + '\n'
                            output_text += 'Change: '+'\t' + ' '+'\n'
                            output_text += 'Measure: ' + '\t' + stab_in.readout(0) + '\n'
                            fout.write(output_text)
                else:
                    ## more than one stabilizers inside the Tableau anti-commute with the input one
                    position = np.where(commutation ==1)[0]
                    position0 = position[-1]                        ## use the last line to make the rest of the stabilizers commute with the input one
                    
                    mtx = self.tableau.tolil()

                    change_list = []
                    for p in position[:-1]:
                        change_list.append(str(p))
                        res = Tableau.tab_add(mtx[position0],mtx[p])
                        if res.count_nonzero() == 0:
                            mtx[p,:]=0
                        else:
                            mtx[p] = res.copy()

                    ## replace the last line by the new stabilizer
                    mtx[position0] = stab_in[0].copy()

                    if not(doc_file is None):                       ## record the action
                        with open(doc_file,'a') as fout:
                            output_text = 'Measure Position: '+'\t' + str(position0) + '\n'
                            output_text = output_text + 'Change: '+'\t' + ' '.join(change_list)
                            output_text += '\t'+self.readout(position0)+'\n'
                            output_text += 'Measure: ' + '\t' + stab_in.readout(0) + '\n'
                            fout.write(output_text)

                    self.tableau = mtx
                    return
        else:
            for j in range(0,stab_in.stab_counts):
                self.measure(stab_in[j],doc_file)
            return

    def permute_qubits(self,indices,doc_file = None):
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
        if len(indices) != self.qubits:
            raise TypeError('The permutation list does not match the qubit number.')

        new_index = np.zeros(2*self.qubits + 1)
        new_index[:self.qubits] = indices
        new_index[self.qubits:-1] = np.array(indices) + self.qubits
        new_index[-1] = 2*self.qubits

        self.tableau = self.tableau.tocoo()
        col_index = np.array([new_index[x] for x in self.tableau.col])
        #print(self.tableau.row)
        #print(row_index)
        self.tableau.col = col_index
        #print(self.tableau.row)
        self.tableau = self.tableau.tocsr()

        if not doc_file is None:
            with open(doc_file,'a') as fout:
                output_text = 'Qubit Permute: \t'
                output_text += ' '.join([str(x) for x in indices])
                output_text += '\n'
                fout.write(output_text)
        return
    
    def _qubit_relabel(self, doc_file = None):
        '''
        This function is to relabel the qubits, to make sure the tableau is in the form of
        [[I  A | B C]
         [0  0 | D E]]
        where I has the maximum dimension.
        '''
        qubit_index = np.zeros(self.qubits)
        x_mtx = self.tableau.tolil()[:,:self.qubits]
        x_mtx = x_mtx.tocsr()
        x_mtx.eliminate_zeros()

        filled_set = set({})
        for j in range(0,self.stab_counts):
            col = x_mtx[j].nonzero()[1]
            # print('j',j)
            if len(col) == 0:
                j -= 1
                break
            else:
                qubit_index[col[0]] = j
                filled_set.add(col[0])
        
        # print('index: ',qubit_index)

        for i in range(0,self.qubits):
            if i in filled_set:
                continue
            else:
                # print('j',j)
                j+=1
                qubit_index[i] = j
        
        # print('index: ',qubit_index)
        self.permute_qubits(qubit_index,doc_file)
        return
    
    def qubit_exchange(self,i,j,doc_file = None):
        '''
        This function exchange the labeling of two qubits
        Input:
            i,j:        The indices for qubit exchange
            doc_file:   The file name for documentation the exchange process
        Output:
            None:       The function will change the self.tableau
        Exception:
            IndexError: If the i, j is invalid, an IndexError will be raised
        '''
        if i>= self.qubits or j>= self.qubits:
            raise IndexError('The input indices are out of range.')
        
        self.tableau = self.tableau.tolil()
        self.tableau[:,i], self.tableau[:,j] = self.tableau[:,j], self.tableau[:,i]     ## exchange X part

        i += self.qubits
        j += self.qubits

        self.tableau[:,i], self.tableau[:,j] = self.tableau[:,j], self.tableau[:,i]     ## exchange Z part

        self.tableau = self.tableau.tocsr()
        return
    
    def _count_nonzero_lines(self, mtx = None):
        '''
        This function count the nonzero lines.
        We assume the tableau has already been reduced.
        Input:
            mtx:            the tableau for the stabilizers
                            If the mtx is None, use self.tableau
        Output:
            line_count:     The number of lines that are not zero (or nontrivial stabilizers)
        '''
        if mtx is None:
            mtx = self.stabilizers()
        mtx = mtx.tocsr()
        mtx.eliminate_zeros()
        row_index = self.stab_counts-1

        while row_index >= 0 and len(mtx[row_index].data) == 0:
            row_index -= 1
        
        return row_index + 1

    def remove_empty_line(self, reduce = False):
        '''
        This function removes the empty lines inside the Tableau
        Input:
            reduce:     whether the tableau needs to be reduced at first, default is False
        Output:
            new_tab:    The new tableau in Tableau class
        
        Note:
            If reduce == False, no reduction process has been done in this function.
            If the self.tableau has not been reduced before, the removal process may not be complete.
        '''

        if reduce:
            self.reduction()
        
        rows = self._count_nonzero_lines()
        new_tab_mtx = self.tableau[:rows]
        new_tab = Tableau(new_tab_mtx)

        return new_tab
    
    def apply_gate(self, gate_name, qubits):
        '''
        A re-written of the gate implementation using the CHP method
        '''
        pass
    
    def gates1(self,gate_name,i):
        '''
        This function perform single-qubit gate on the qubit.
        Input:
            gate_name:  The gate name for the single qubit gate
                        gate name should be: [X, Y, Z, H, S] (case insensitive)
            i:          The index for the qubit that needs to perform the single-qubit gate
        Output:
            None:       This will change the self.tableau
        Exception:
            ValueError: If the gate name is invalid, a ValueError will be raised.
            IndexError: If the qubit number is out of the possible qubit number
        '''
        if not gate_name in SINGLE_GATE_SET:
            raise ValueError('The gate operation is invalide.')
            
        if i >= self.qubits:
            raise IndexError('The qubit index is out of range.')

        mtx = self.tableau.tolil()
        if gate_name == 'x' or gate_name == 'X':
            ## Pauli X gate:
            ##      X -> X
            ##      Z -> -Z
            ##      Y -> -Y
            _ , row_index = self.search_nonzero_rows(i+self.qubits)     ## row_index is the ones in the Z segment.
            for j in row_index:
                mtx[j,-1] = ~mtx[j,-1]      ## reverse the sign of these rows
            self.tableau = mtx
            return
        elif gate_name == 'z' or gate_name == 'Z':
            ## Pauli Z gate:
            ##      X -> -X
            ##      Z -> Z
            ##      Y -> -Y
            _ , row_index = self.search_nonzero_rows(i)     ## row_index is the ones in the X segment.
            for j in row_index:
                mtx[j,-1] = ~mtx[j,-1]      ## reverse the sign of these rows
            self.tableau = mtx
            return
        elif gate_name == 'y' or gate_name == 'Y':
            ## Pauli Y gate:
            ##      X -> -X
            ##      Z -> -Z
            ##      Y -> Y
            _ , row_index_x = self.search_nonzero_rows(i)               ## row_index is the ones in the X segment.
            _ , row_index_z = self.search_nonzero_rows(i+self.qubits)   ## row_index is the ones in the Z segment.
            row_index_x = set(row_index_x)
            row_index_z = set(row_index_z)
            
            ## the following row_index_xz is for the row indices that are either X or Z, not Y.
            ## the way we do it is: 
            row_index_xz = row_index_x.symmetric_difference(row_index_z)

            for j in row_index_xz:
                mtx[j,-1] = ~mtx[j,-1]      ## reverse the sign of these rows, (these rows have X or Z, not Y)
            self.tableau = mtx
            return
        elif gate_name == 'h' or gate_name == 'H':
            ## Hadamard gate:
            ##      X -> Z
            ##      Z -> X
            ##      Y -> -Y
            _ , row_index_x = self.search_nonzero_rows(i)               ## row_index is the ones in the X segment.
            _ , row_index_z = self.search_nonzero_rows(i+self.qubits)   ## row_index is the ones in the Z segment.
            row_index_x = set(row_index_x)
            row_index_z = set(row_index_z)

            ## the following row_index_y is for the row indices that are Y.
            ## the way we do it is: 
            row_index_y = row_index_x.intersection(row_index_z)

            ## exchange X and Z columns, to achive X <->Z, Y -> Y
            mtx[:,i], mtx[:, i+self.qubits] = mtx[:, i+self.qubits], mtx[:, i]

            ## reverse the sign of Y rows
            for j in row_index_y:
                mtx[j,-1] = ~mtx[j,-1]          ## reverse sign
            self.tableau = mtx
            return
        elif gate_name == 's' or gate_name == 'S':
            ## phase gate:
            ##      X -> Y
            ##      Y -> -X
            ##      Z -> Z
            _ , row_index_x = self.search_nonzero_rows(i)               ## row_index is the ones in the X segment.

            for j in row_index_x:
                if mtx[j,i+self.qubits] == 1:
                    ## means current stabilizer is Y, not X
                    ## we need to reverse the sign bit for this stabilizer
                    mtx[j,-1] = ~mtx[j,-1]
                mtx[j,i+self.qubits] = ~mtx[j,i+self.qubits]        ## exchange X <-> Y for this qubit
            self.tableau = mtx
            return
    
    def standard_form(self,doc_file = None):
        '''
        This function reduce the tableau to the standard form.
        This function should be useful for the stabilizers of error correction code.
        For error correction code, the stabilizer count < qubit number.
            The standard form of the stabilizer code is
            [[I A1 A2 | B 0 C]
             [0 0  0  | D I E]]
        If the stabilizer is full rank (stabilizer count = qubit number), the stabilizer tableau can be reduced to
            [[I A | B 0]
             [0 0 | D I]]
        '''

        ## the following code make sure that the stabilizer tableau is in the form of
        ## [[I A | B C]
        ##  [0 0 | D E]]
        self.reduction(doc_file)

        ##########
        # Till here, the X part has already been transformed to have the largest identity block matrix.
        # Next, we want to reduce the Z part below.
        ##########

        ## extract the size of identity block in X part
        x_rank = self._count_nonzero_lines(self.tableau[:,:self.qubits])

        if x_rank == self.stab_counts:
            ## if the x_rank is equal to the length of the stabilizers, then we are done.
            return

        row_pointer = x_rank      ## current row position is the first zero row for X
        column_pointer = self.qubits + x_rank -1   ## current column position, the first column in the desired I in Z part

        for column_pointer in range(self.qubits + x_rank,self.qubits*2):
            #print(self.tableau.todense())
            #print(row_pointer)
            _,ones = self.search_nonzero_rows(column_pointer,start_row = row_pointer)
            ## group the rest of the rows by either the column == 1 or not
            if len(ones) == 0:
                ## there is no rows that has "1" on "column_pointer" column position
                ## we skip this column and do nothing
                continue
            else: 
                #print('reduce down ones:',ones)
                self._reduce_down(ones, doc_file)    ## make sure the row: ones[0] is the only row that have a 1 on the "column_pointer" column position
                self.row_insert(ones[0],row_pointer,doc_file)  ## insert this row to the current row_pointer pointed rows

                row_pointer += 1
                
                #print(self.tableau.todense())
                
                ## following: we try to use this new row, to elimitate the previous rows that has one on this working column
                _, ones_up = self.search_nonzero_rows(column_pointer,final_row = row_pointer)
                self._reduce_up(ones_up,doc_file)
                
                #print(self.tableau.todense())
                
                self.tableau.eliminate_zeros()
                
                ## if we are working on the last row, we are done
                if row_pointer == self.tableau.shape[0]:
                    break

        ## Here we want to exchange qubits in the rest of the qubits to make the I in Z part of the tableau to be the maximum size
        qubit_index = np.zeros(self.qubits)
        qubit_index[:x_rank] = range(x_rank)

        ## Here we only focus on the [I, E] part in the Z-part of the Tableau
        z_mtx = self.tableau.tolil()[:,self.qubits+x_rank:-1]
        z_mtx = z_mtx.tocsr()
        z_mtx.eliminate_zeros()

        filled_set = set(range(x_rank))
        for j in range(x_rank,self.stab_counts):
            col = z_mtx[j].nonzero()[1]
            #print('col',col, 'j', j)
            if len(col) == 0:
                j -= 1
                break
            else:
                qubit_index[col[0]+x_rank] = j
                filled_set.add(col[0]+x_rank)
        
        # print(qubit_index)

        for i in range(x_rank,self.qubits):
            if i in filled_set:
                continue
            else:
                j+=1
                qubit_index[i] = j
        
        # print(qubit_index)
        
        self.permute_qubits(qubit_index,doc_file)
        return
    
    def turn_to_graph(self,doc_file=None, gate_record=False):
        '''
        This function turn the stabilizer state defined by the Tableau to the equivalent graph states
        If the Tableau is not full-rank (rank != qubit number):
        the function will fail.
        Input:
            This function uses the current instance's tableau
            doc_file:       The document for recording the operations
            gate_record:    if the gate applied to turn to the graph state is returned or not.
                            Default is False
        Output:
            success:       True/False
            
            If success, the graph state tableau will be returned as a new Tableau class instance
            new_graph:      The new stabilizer state tableau
        
        Notice that:
            This function will automatically change the current tableau to the standard form.
            current version (2020/12/16):
                Ignore the extra overall phase

            If the stabilizer is full rank (stabilizer count = qubit number), the stabilizer tableau can be reduced to
            [[I A | B 0]
             [0 0 | D I]]

            If B has any nonzero diagonal elements, apply S gates to turn Y -> X
            ==> 
            [[I A | B' 0]
             [0 0 | D  I]]
            where B' has zero diagonal elements
            
            ==> Apply Hadamard on the second block of qubits to exchange X and Z
            [[I 0 | B' A]
             [0 I | D  0]]
            
            This is a standard representation of the graph state, upto overall phase on some stabilizers (rows)

            ==> Apply Z to those stabilizers with a (-) sign to correct the sign

            This function, we apply a short-cut, to virtually apply the gates.
        '''
        graph_tableau = Tableau(self.tableau)
        graph_tableau.standard_form(doc_file)

        if graph_tableau._count_nonzero_lines() !=graph_tableau.qubits:
            return False
        
        ## extract the size of identity block in X part
        x_rank =graph_tableau._count_nonzero_lines(graph_tableau.tableau[:,:graph_tableau.qubits])

        graph_mtx = graph_tableau.tableau
        #print(graph_mtx.todense())

        graph_mtx = graph_mtx.tolil()
        i1 = graph_mtx[:x_rank,:x_rank]
        a_mtx = graph_mtx[:x_rank,x_rank:self.qubits]
        b_mtx = graph_mtx[:x_rank,self.qubits:self.qubits+x_rank]
        sign1 = sparse.csr_matrix(np.zeros([x_rank,1]))

        d_mtx = graph_mtx[x_rank:,self.qubits:self.qubits+x_rank]
        i2 = graph_mtx[x_rank:,self.qubits+x_rank:-1]
        sign2 = sparse.csr_matrix(np.zeros([self.qubits-x_rank,1]))

        ## set the diagonal elements of b_mtx to be zero
        s_gate = []
        for i in range(b_mtx.shape[0]):
            # record the s_gate applied on the qubit
            if b_mtx[i,i] != 0:
                s_gate.append(i)
            b_mtx[i,i] = 0
        
        if gate_record:
            i2_dim = i2.shape[0]
            h_gate = [self.qubits - (i+1) for i in range(i2_dim)]

            z_gate = []
            for i in range(graph_mtx.shape[0]):
                if graph_mtx[i,-1] != 0:
                    z_gate.append(i)
        
        ##############################
        ##
        ## debug code:
        ## print("a_mtx dim: ",a_mtx.shape)
        ## print("b_mtx dim: ",b_mtx.shape)
        ## print("d_mtx dim: ",d_mtx.shape)
        ## print("i1 dim: ",i1.shape)
        ## print("i2 dim: ",i2.shape)
        ## print("sign1 dim: ",sign1.shape)
        ## print("sign2 dim: ",sign2.shape)
        ##############################
        
        graph_mtx_new = sparse.bmat([[i1,None,b_mtx,a_mtx,sign1],[None,i2,d_mtx,None,sign2]])
        graph_mtx_new = graph_mtx_new.tocsr()
        graph_mtx_new.eliminate_zeros()
        graph_tableau = Tableau(graph_mtx_new)
        if gate_record:
            return graph_tableau,s_gate,h_gate,z_gate
        else:
            return graph_tableau

    def entanglement_entropy(self,subsystem):
        '''
        Solving the entanglement entropy of the subsystem
        Input:
            subsystem:  a list (or ndarray) which saves the indices of the qubits that belong to the sub-system we are going to calculate.
        Output:
            entropy:    The Von Newmann entropy for the subsystem.
        '''
        ## Step 1: permute the qubit labeling to gather all the qubits inside the sub-system to the right.
        qubit_index = [-1] * self.qubits
        pos = self.qubits - 1
        for i in range(len(subsystem),0,-1):
            qubit_index[subsystem[i-1]] = pos
            pos -= 1
        
        for i in range(len(qubit_index),0,-1):
            if qubit_index[i-1] == -1:
                qubit_index[i-1] = pos
                pos -= 1
        # print('qubit_index = ',qubit_index)
        temp_tableau = Tableau(self.tableau)
        temp_tableau.permute_qubits(qubit_index)

        ## Here the phase bit is not important.
        ## So we just extract the stabilizers as a binary matrix
        stab_mtx =temp_tableau.stabilizers()
        stab_mtx = BM(stab_mtx)


        ##########
        # print('stab_mtx = ')
        # print(stab_mtx)
        ##########
        ## Step 2: reformulate the tableau matrix by making the columns as [x1, z1, x2, z2 ...]
        new_col_index = list(range(2*temp_tableau.qubits))
        new_col_index[:temp_tableau.qubits] = list(range(0,2*temp_tableau.qubits,2))
        new_col_index[temp_tableau.qubits:] = list(range(1,2*temp_tableau.qubits,2))

        ## Step 3: Gaussian reduction
        stab_mtx.col_permute(new_col_index)
        stab_mtx.reduction()

        ##########
        # print('permute: ',new_col_index)
        # print('stab_mtx: ')
        # print(stab_mtx)
        ##########

        ## Step 4: extract the stabilizers that only acts on the subsystem qubits
        # The sub-system position is at 2*[self.qubits - len(sybsystem)]
        # so we want the rows that pivot elements is the sub-system columns.
        r_index = -1
        for row_num,row in enumerate(stab_mtx.bm):
            col = row.nonzero()[1]
            col.sort()
            if len(col) != 0:
                if col[0] >= 2*(temp_tableau.qubits - len(subsystem)):
                    r_index = row_num
                    break
        
        ##########
        # print('qubit num = ',len(subsystem))
        # print('r_index = ',r_index)
        ##########

        if r_index == -1:
            return len(subsystem)
        else:
            temp_mtx =stab_mtx.bm[r_index:,2*(temp_tableau.qubits - len(subsystem)):]
            temp_mtx = BM(temp_mtx) 
            r_sub = temp_mtx.rank(False)
            #print('r_sub = ',r_sub)
            return len(subsystem) - r_sub

    def height(self, sep_index, reorder=None):
        '''
        solving the height function of the stabilizer state when bi-partition the system at qubit sep_index.
        The bi-partition of the system is [0,1, ..., sep_index-1] vs [sep_index, ... n]
        Input:
            sep_index:      the seperation index of the qubit to perform bi-partition of the system
            reorder:        If this is given, we will relabel the qubit as given.
                            e.g., if the reorder list is given as [0,2,3,1] for a four-qubit system
                            the new qubit 1 is the old qubit 2, new qubit 2 is the old qubit 3, etc.
        Output:
            h:              the height function value with the bipartition
        
        Note that when reorder is given, the self instance is re-ordered. 
        '''

        ## step 1: permute the qubits with the given order, if there is any.
        if reorder is not None:
            self.permute_qubits(reorder)
        
        ## step 2: get the stabilizers
        ## Note that the signs of the stabilizers does not affect the final calculation
        ## we will use the binary matrix instead.

        stab_temp = self.stabilizers()
        stab_temp = stab_temp.toarray()
        bm = BM(stab_temp)
        ## permute the columns in the binary matrix to make X and Z column for a single qubit togather.
        index_list = [2*i for i in range(self.qubits)] + [2*i+1 for i in range(self.qubits)]
        bm.col_permute(index_list)
        bm.reduction()

        ## this is the part of the sub-system which is to the left of the sep_index qubit
        bm_chop1 = bm[:, :2*sep_index]
        stab_num = bm_chop1.non_zero_lines()

        ## this is the part of the RREF part of the tableau
        bm_chop = bm[stab_num:,2*sep_index:]
    
        rank = bm_chop.rank(False)
        #print(stab_num, rank)
        return self.qubits - sep_index - rank

    def height_sweep(self, reorder=None):
        '''
        solving the height function of the stabilizer as we sweep the bipartition site in the given system.
        if reorder is given, the qubit index will be relabled according to the list given in reorder.
        Input:
            reorder:        relabel the qubits as given.
        Output:
            h_list:         the list of height values with different bi-partition
        '''

        ## step 1: permute the qubits with the given order, if there is any.
        if reorder is not None:
            self.permute_qubits(reorder)

        ## step 2: get the stabilizers

        stab_temp = self.stabilizers()
        stab_temp = stab_temp.toarray()
        bm = BM(stab_temp)
        ## permute columns in the binary matrix
        index_list = [2*i for i in range(self.qubits)] + [2*i+1 for i in range(self.qubits)]
        bm.col_permute(index_list)
        bm.reduction()

        ## step 3: sweep the bi-partition positions
        h_list = []
        for p in range(0,self.qubits+1):
            bm_chop1 = bm[:, :2*p]
            stab_num = bm_chop1.non_zero_lines()

            bm_chop = bm[stab_num:,2*p:]
    
            rank = bm_chop.rank(False)
            h_list.append(self.qubits - p - rank)
        
        return h_list



def graph_equivalence(adj_mtx1, adj_mtx2):
    '''
    This function uses the two adjacency matrices to determine if the two graphs are LC equivalent.
    The method is from:
    (1) M. Van den Nest et al, Phys. Rev. A, 69, 022316 (2004)
    (2) M. Van den Nest et al, Phys. Rev. A, 70, 034302 (2004)
    Input:
        adj_mtx1, adj_mtx2:     Two adjacency matrices for two graphs
    Output:
        equivalence:        True/False

        If the two graphs are equivalent, the corresponding transformation matrix will be given.
        q_mtx:              The corresponding transformation matrix to transform these two graphs to be the same.
    '''
    ## convert the two adjacency matrices to BinaryMatrix class instances
    adj1 = BM(adj_mtx1)
    adj2 = BM(adj_mtx2)

    ## for valid adjacency matrices, the two binary matrices should be square matrices
    if len(adj1.shape) != 2:
        raise IndexError('The adj_mtx1 had wrong dimension as a graph adjacency matrix.')
    elif adj1.shape[0] != adj1.shape[1]:
        raise IndexError('The adj_mtx1 had wrong dimension as a graph adjacency matrix.')

    if len(adj2.shape) != 2:
        raise IndexError('The adj_mtx2 had wrong dimension as a graph adjacency matrix.')
    elif adj2.shape[0] != adj2.shape[1]:
        raise IndexError('The adj_mtx2 had wrong dimension as a graph adjacency matrix.')
    
    ## if two graphs have different vertices number, the two graphs are definitely in-equivalent.
    if adj1.shape[0] != adj2.shape[0]:
        return False
    
    ## if the two adjacency matrices are the same, the two graphs are equivalent.
    test = adj1 + adj2
    test.bm.eliminate_zeros()
    if len(test.bm.data) == 0:
        return (True,sparse.eye(4*adj1.shape[0]))
    
    ## if the easy cases cannot tell, we will need to solve the transformation matrices
    ## First step:
    ##      Construct the linear equations to make it ready for solving in GF(2)
 
    ## the coefficient matrix is a (n^2 * 4n) matrix
    ## the variables are aligned as {a's, b's, c's, d's} with each sub-vector is length n.
    vn = adj1.shape[0] ## vertices number
    elem = []
    row = []
    col = []
    for j in range(0,vn):
        for k in range(0,vn):
            row_num = j*vn + k 

            ## a's coefficients
            elem.append(adj1.bm[j,k])
            row.append(row_num)
            col.append(k)

            ## b's coefficients
            if j == k:
                elem.append(1)
                row.append(row_num)
                col.append(vn+k)
            
            ## d's coefficients
            elem.append(adj2.bm[j,k])
            row.append(row_num)
            col.append(3*vn+j)

            ## c's coefficients
            for m in range(0,vn):
                #elem.append(adj1.bm[j,m] * adj2.bm[m,k])
                elem.append(adj1.bm[m,j] * adj2.bm[m,k])
                row.append(row_num)
                col.append(2*vn+m)
    
    co_mtx = sparse.csr_matrix((elem,(row,col)),shape=(vn**2,4*vn))
    co_mtx = co_mtx.astype(int)

    co_BM = BM(co_mtx)

    ## use the BinaryMatrix reduction method to reduce the coefficient matrix (co_BM) to RREF and solve the null-space
    null_space = co_BM.nullspace()

    ## if the co-efficient matrix has rank = 4n, i.e, the null-space contains a single zero vector,
    ## this means the two graphs are not equivalent.
    null_space.bm.eliminate_zeros()
    if len(null_space.bm.data) == 0:
        return False
    
    ## if there are non-trivial vectors inside the null-space
    ## we will test the constrain AD + BC == 1
    ## using the lamma in Ref(2) as Bikun's code did. 

    ## test if the Q matrix has full rank in GF(2) for each null-space vectors.
    ## if the dimension of the null-space is lower than 4
    ## we will test all possible combinations of valid Q matrices
    
    ## if the dimension of null-space is larger than 4
    ## we try to select two valid Q matrices, and add them up, to test their rank

    final_q_mtx = None

    null_dim = null_space.shape[0]
    if null_dim <= 4:
        ## case 1
        for i in range(0,2**null_dim):
            ## enumerate all possible combinations

            num = i
            test_vec = BM(np.zeros([1,null_space.shape[1]]).astype(int))
            for j in range(0,null_dim):
                if num%2 == 1:
                    test_vec = test_vec + null_space[j]
                num = num >> 1
            
            q_mtx = _get_Q_mtx(test_vec,vn)
            q_mtx_bm = BM(q_mtx)
            if q_mtx_bm.rank(False) == 2*vn:
                final_q_mtx = q_mtx
                break
    
        if not final_q_mtx is None:
            return (True,final_q_mtx)
        else:
            return False
    else:
        ## case 2, only test upto 2-sum of null-vectors
        ## test single null-vector
        for i in range(0,null_dim):
            q_mtx = _get_Q_mtx(null_space[i],vn)
            q_mtx_bm = BM(q_mtx)
            if q_mtx_bm.rank(False) == 2*vn:
                return (True,q_mtx)
            
        for i in range(0,null_dim):
            for j in range(0,null_dim):
                test_vec = null_space[i] + null_space[j]
                q_mtx = _get_Q_mtx(test_vec,vn)
                q_mtx_bm = BM(q_mtx)
                if q_mtx_bm.rank(False) == 2*vn:
                    final_q_mtx = q_mtx
                    break
            if not (final_q_mtx is None):
                break
        
        if final_q_mtx is None:
            return False
        else:
            return (True,final_q_mtx)

def _get_Q_mtx(null_vector,vn = None):
    length = null_vector.shape[1]
    if isinstance(null_vector,BM):
        null_vector = null_vector.bm
    
    if vn is None:
        vn = length//4

    mtx = null_vector.toarray()
    a_diag = mtx[0,:vn]
    b_diag = mtx[0,vn:2*vn]
    c_diag = mtx[0,2*vn:3*vn]
    d_diag = mtx[0,3*vn:]

    a_mtx = sparse.diags(a_diag)
    b_mtx = sparse.diags(b_diag)
    c_mtx = sparse.diags(c_diag)
    d_mtx = sparse.diags(d_diag)

    q_mtx = sparse.bmat([[d_mtx,c_mtx],[b_mtx,a_mtx]])
    q_mtx = q_mtx.astype(int)
    return q_mtx
    