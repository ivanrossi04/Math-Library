template <typename T>
Matrix<T>::Matrix(){
    this -> data = nullptr;
    this -> rows = 0;
    this -> columns = 0;
}

template <typename T>
Matrix<T>::Matrix(size_t n, size_t m, T val){
    this -> rows = n;
    this -> columns = m;

    // the matrix gets filled cell by cell with a given value
    // dim is the total number off cells in the matrix
    size_t dim = n * m;
    this -> data = new T[dim];
    for(size_t i = 0; i < dim; i++){
            this -> data[i] = val;
    }
}

template <typename T>
Matrix<T>::Matrix(const std::initializer_list<std::initializer_list<T>>& matrix){
    this -> rows = matrix.size();

    if(this -> rows){
        this -> columns = matrix.begin() -> size();
        
        this -> data = new T[this -> rows * this -> columns];

        size_t i = 0;
        for(std::initializer_list<T> row: matrix){
            for(T element: row){
                this -> data[i] = element;
                i++;
            }
        }
    } else{
        this -> columns = 0;
        this -> data = nullptr;
    }
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& m){
    this -> rows = m.rows;
    this -> columns = m.columns;

    size_t dim = m.rows * m.columns;
    this -> data = new T[dim];

    for(size_t i = 0; i < dim; i++)
        this -> data[i] = m.data[i];
    
}

template <typename T>
void Matrix<T>::operator=(const Matrix<T>& m){
    // before a new assignment it is necessary to check for old data to delete
    if(this -> data != nullptr) {
        delete[] this -> data;
        this -> data = nullptr;
    }

    *this = Matrix<T>(m);
}

template <typename T>
void Matrix<T>::operator=(const std::initializer_list<std::initializer_list<T>>& matrix){
    // before a new assignment it is necessary to check for old data to delete
    if(this -> data != nullptr) {
        delete[] this -> data;
        this -> data = nullptr;
    }

    *this = Matrix<T>(matrix);
}

template <typename T>
void Matrix<T>::operator=(Matrix<T>&& m){
    this -> rows = m.rows;
    this -> columns = m.columns;

    if(this -> data != nullptr) {
        delete[] this -> data;
        this -> data = nullptr;
    }

    this -> data = m.data;

    m.data = nullptr;
}

// ----------------------------------------------------------------------
template<typename T>
Matrix<T> Matrix<T>::identity(size_t n){
    // matrix gets declared
    Matrix<T> id;

    id.rows = n;
    id.columns = n;

    id.data = new T[n * n];

    // double iteration is needed to easily check for a diagonal element
    for(size_t i = 0; i < n; i++)
        for(size_t j = 0; j < n; j++)
            id.data[i * id.columns + j] = (T)(i == j); 
    

    return id;
}

template<typename T>
T* Matrix<T>::operator[](size_t index){
    if(index > this -> rows) throw "Index out of bounds";
    return &(this -> data[index * this -> columns]);
}

template<typename T>
size_t Matrix<T>::getRows(){return this -> rows;}

template<typename T>
size_t Matrix<T>::getColumns(){return this -> columns;}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& m){
    if(this -> rows != m.rows || this -> columns != m.columns) throw "Invalid operation: can't sum matrices of different dimensions";

    Matrix result(this -> rows, this -> columns);
    size_t dim = result.rows * result.columns;

    // addition is performed between the two matrices cell by cell
    for(size_t i = 0; i < dim; i++)
        result.data[i] = this -> data[i] + m.data[i];

    return result;
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& m){
    if(this -> rows != m.rows || this -> columns != m.columns) throw "Invalid operation: can't sum matrices of different dimensions";

    size_t dim = this -> rows * this -> columns;

    for(size_t i = 0; i < dim; i++)
        this -> data[i] += m.data[i];

    return *(this);
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& m){
    if(this -> rows != m.rows || this -> columns != m.columns) throw "Invalid operation: can't sum matrices of different dimensions";

    Matrix result(this -> rows, this -> columns);
    size_t dim = result.rows * result.columns;

    // subtraction is performed between the two matrices cell by cell
    for(size_t i = 0; i < dim; i++)
        result.data[i] = this -> data[i] + m.data[i];

    return result;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& m){
    if(this -> rows != m.rows || this -> columns != m.columns) throw "Invalid operation: can't sum matrices of different dimensions";

    size_t dim = this -> rows * this -> columns;

    for(size_t i = 0; i < dim; i++)
        this -> data[i] -= m.data[i];

    return *(this);
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& m){
    if(this -> columns != m.rows) throw "Invalid operation: operation between matrices (A*B) must follow the rule of dimensions A: m*n, B:n*p";
    
    /*
    given two matrices A(n x m) and B(m x p) 
    the multiplication result C of A x B is a matrices with n rows and p columns
    where the C_ij element is given by: C_ij = sum{Aik * Bkj}
    */

    Matrix result = Matrix(this -> rows, m.columns);

    // iteration on the n rows of the first matrix
    for(size_t i = 0; i < this -> rows; i++){
        // iteration on the p columns of the second matrix
        for(size_t j = 0; j < m.columns; j++){
            // computation of the result value in the C_ij cell
            for(size_t k = 0; k < this -> columns; k++){
                result.data[i * m.columns + j] += this -> data[i * this -> columns + k] * m.data[k * m.columns + j];
            }
        }
    }

    return result;
}

template<typename T>
Matrix<T> operator*(const T& val, const Matrix<T>& m){
    Matrix<T> result(m);

    for(size_t i = 0; i < result.rows * result.columns; i++)
        result.data[i] *= val;

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T& val){

    Matrix<T> result(*(this));

    for(size_t i = 0; i < result.rows * result.columns; i++)
        result.data[i] *= val;

    return result;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const T& val){
    for(size_t i = 0; i < this -> rows * this -> columns; i++)
        this -> data[i] *= val;

    return *(this);
}

template<typename T>
Matrix<T> Matrix<T>::transpose(const Matrix<T>& m){
    Matrix transposed = Matrix(m.columns, m.rows);

    for(size_t i = 0; i < m.rows; i++){
        for(size_t j = 0; j < m.columns; j++){
            transposed.data[i * m.rows + j] = m.data[j * m.columns + i];
        }
    }

    return transposed;
}

template<typename T>
Matrix<T> Matrix<T>::invert(Matrix<T> m){
    Matrix<T> inverse = Matrix<T>::identity(m.columns);

    // every operation on the matrix m is repeated on the inverse matrix

    for(size_t i = 0; i < m.rows; i++){
        if(m.data[i * m.columns + i] == 0){
            size_t j = i + 1;

            // search pivot on the current column
            bool found_pivot = false;
            while(!found_pivot && j < m.columns){
                if(m.data[j * m.columns + i] != 0){
                    inverse.swapRows(i, j);

                    m.swapRows(i, j);
                    found_pivot = true;
                }
                j++;
            }

            // if a pivot is not found the matrix is singular
            if(!found_pivot) throw "Error: The given matrix is singular (non invertible)";
        }

        // operations on the rows below the current to eliminate the cells below the pivot
        for(size_t j = 0; j < m.rows; j++){
            if(j != i){
                inverse.addRow(j, i, -(m.data[j * m.columns + i] / m.data[i * m.columns + i]));
                m.addRow(j, i, -(m.data[j * m.columns + i] / m.data[i * m.columns + i]));
            }
        }

        inverse.multiplyRow(i, 1.0 / m.data[i * m.columns + i]);
        m.multiplyRow(i, 1.0 / m.data[i * m.columns + i]);
    }

    return inverse;
}

template<typename T>
T Matrix<T>::det(const Matrix<T>& m){
    if(m.rows != m.columns) throw "Invalid operation: can't compute the determinant of a non square matrix";
    else if(m.rows == 0) return 0;
    else if(m.rows == 1) return m.data[0];
    else if(m.rows == 2) return (m.data[0] * m.data[3]) - (m.data[1] * m.data[2]);
    else if(m.rows == 3) {
        // static implementation of the laplacian formula
        return  m.data[0] * (m.data[4] * m.data[8] - m.data[7] * m.data[5]) -
                m.data[1] * (m.data[3] * m.data[8] - m.data[6] * m.data[5]) +
                m.data[2] * (m.data[3] * m.data[7] - m.data[6] * m.data[4]);
    } else{
        // implementation of gaussian elimination
        Matrix<T> copy(m);
        T det = 1;

        for(size_t i = 0; i < copy.rows - 1; i++){
            if(copy.data[i * copy.columns + i] == 0){
                size_t j = i + 1;

                // search pivot on the current column
                bool found_pivot = false;
                while(!found_pivot && j < copy.columns){
                    if(copy.data[j * copy.columns + i] != 0){
                        copy.swapRows(i, j);
                        found_pivot = true;
                        det *= -1;
                    }
                    j++;
                }

                // if a pivot is not found the matrix is singular
                if(!found_pivot) return 0;
            }

            // the pivot value on the diagonal is multiplied
            det *= copy.data[i * copy.columns + i];

            // operations on the rows below the current to eliminate the cells below the pivot
            for(size_t j = i + 1; j < copy.rows; j++){
                copy.addRow(j, i, -(copy.data[j * copy.columns + i] / copy.data[i * copy.columns + i]));
            }
        }
        return det * copy.data[m.rows * m.rows - 1];
    }
}

template <typename U>
std::ostream& operator <<(std::ostream &out, const Matrix<U> &m) {

    /* print example: std::cout << Matrix({{2, 1},{0, -2}});
    2   1
    0   -2

    */ 
    for(size_t i = 0; i < m.rows; i++){
        for(size_t j = 0; j < m.columns; j++) out << (m.data[i * m.columns + j]) << ((j != m.columns - 1) ? "\t" : "");
        out << "\n";
    }

    return out;
}

template <typename T>
Matrix<T>::~Matrix(){
    if(this -> data != nullptr){
        delete[] this -> data;
        this -> data = nullptr;
    }
}

// ------- Linear algebra functions implementation -------
template<typename T>
void Matrix<T>::swapRows(size_t row1, size_t row2){

    if(row1 >= this -> rows || row2 >= this-> rows) throw "Index out of bounds";

    T temp;
    for(size_t i = 0; i < this -> columns; i++) {
        temp = this -> data[row1 * this -> columns + i];
        this -> data[row1 * this -> columns + i] = this -> data[row2 * this -> columns + i];
        this -> data[row2 * this -> columns + i] = temp;
    }
}

template<typename T>
void Matrix<T>::multiplyRow(size_t row, T val){
    if(row >= this -> rows) throw "Index out of bounds";

    for(size_t i = 0; i < this -> columns; i++)
        this -> data[row * this -> columns + i] *= val;
}

template<typename T>
void Matrix<T>::addRow(size_t row1, size_t row2, T val){
    if(row1 >= this -> rows || row2 >= this -> rows) throw "Index out of bounds";

    for(size_t i = 0; i < this -> columns; i++)
        this -> data[row1 * this -> columns + i] += val * this -> data[row2 * this -> columns + i];
}

template<typename T>
Matrix<T>* LPU_decomposition(const Matrix<T> &m){

    Matrix<T>* LPU = new Matrix<T>[3];

    // initialize U as a copy of m
    LPU[2] = m;


    size_t n = LPU[2].getRows();
    if(n != LPU[2].getColumns()) throw "Invalid format: need a square matrix to perform operation";

    // initialize an L an P as identity matrices
    Matrix<T> Id = Matrix<T>::identity(n);
    LPU[0] = Id;
    LPU[1] = Id;

    for(size_t i = 0; i < n - 1; i++){
        if(LPU[2][i][i] == 0){
            size_t j = i + 1;
            bool found_pivot = false;
            while(!found_pivot && j < n){
                // search pivot on the current column

                if(LPU[2][i][j] != 0){
                    // when a pivot is found:

                    // the row with the pivot gets swapped with the current row in the matrix U
                    LPU[2].swapRows(i, j);

                    // the corresponding rows on P get swapped too
                    LPU[1].swapRows(i, j);

                    // in the matrix L the non diagonal cells of the corresponding rows get swapped
                    LPU[0] -= Id;
                    LPU[0].swapRows(i, j);
                    LPU[0] += Id;

                    found_pivot = true;
                }
                j++;
            }

            if(!found_pivot) throw "Error: matrix isn't invertible";
        }

        for(size_t j = i + 1; j < n; j++){
            // the ratio value between rows is assigned to the corresponding cell
            LPU[0][j][i] = LPU[2][j][i] / LPU[2][i][i];

            // the j row gets eliminated
            LPU[2].addRow(j, i, -LPU[0][j][i]);
        }
    }

    return LPU;
}

template<typename T>
Matrix<T>* QR_decomposition(const Matrix<T> &m){
    Matrix<T>* QR = new Matrix<T>[2];

    // Q matrix is initialized as the transposed of m (this simplifies the vector declaration)
    QR[0] = Matrix<T>::transpose(m);

    size_t n = QR[0].getRows();
    if(n != QR[0].getColumns()) throw "Invalid format: need a square matrix to perform operation";   

    // the algorithm uses the gram-schmidt procedure
    Vec<T> *u = new Vec<T>[n]; 
    for(size_t i = 0; i < n; i++){

        // the current row gets redeclared as a vector
        u[i] = Vec<T>(QR[0][i], n);
        
        // the current operating vector gets orthogonalized in respect to the previous ones
        for(size_t j = 0; j < i; j++)
            u[i] -= Vec<T>::dot(u[i], u[j]) * u[j];

        // the current vector gets normalized
        u[i] *= sqrt(1 / (double) u[i].getNorm());

        // the current vector (orthonormalized) is stored as the new row for Q(transposed)
        for(size_t j = 0; j < n; j++) QR[0][i][j] = u[i][j];
    }
    
    QR[1] = QR[0] * m;
    QR[0] = Matrix<T>::transpose(QR[0]);

    return QR; 
}