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

    this -> data = new T*[n];
    for(size_t i = 0; i < n; i++){
        this -> data[i] = new T[m];
        for(size_t j = 0; j < m; j++)
            this -> data[i][j] = val;
    }
}

template <typename T>
Matrix<T>::Matrix(const std::initializer_list<std::initializer_list<T>>& matrix){
    this -> rows = matrix.size();

    if(this -> rows){
        this -> columns = matrix.begin() -> size();
        
        this -> data = new T*[this -> rows];

        size_t i = 0;
        for(std::initializer_list row: matrix){
            this -> data[i] = new T[this -> columns];
            std::copy(row.begin(), row.end(), this -> data[i]);
            i++;
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
    this -> data = new T*[this -> rows];

    for(size_t i = 0; i < this -> rows; i++){
        this -> data[i] = new T[this -> columns];
        for(size_t j = 0; j < this -> columns; j++){
            this -> data[i][j] = m.data[i][j];
        }
    }
}

template <typename T>
void Matrix<T>::operator=(const Matrix<T>& m){
    if(this -> data != nullptr) {
        for(size_t i = 0; i < this -> rows; i++) {
            delete[] this -> data[i];
            this -> data[i] = nullptr;
        }
        delete[] this -> data;
        this -> data = nullptr;
    }

    this -> rows = m.rows;
    this -> columns = m.columns;
    this -> data = new T*[this -> rows];
    for(size_t i = 0; i < this -> rows; i++){
            this -> data[i] = new T[this -> columns];
            for(size_t j = 0; j < this -> columns; j++){
                this -> data[i][j] = m.data[i][j];
            }
        }
}

template <typename T>
void Matrix<T>::operator=(const std::initializer_list<std::initializer_list<T>>& matrix){
    if(this -> data != nullptr) {
        for(size_t i = 0; i < this -> rows; i++) {
            delete[] this -> data[i];
            this -> data[i] = nullptr;
        }
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
        for(size_t i = 0; i < this -> rows; i++) {
            delete[] this -> data[i];
            this -> data[i] = nullptr;
        }
        delete[] this -> data;
    }

    this -> data = m.data;

    m.data = nullptr;
}

template<typename T>
Matrix<T> Matrix<T>::identity(size_t n){

    Matrix<T> id;

    id.rows = n;
    id.columns = n;

    id.data = new T*[n];
    for(size_t i = 0; i < n; i++){
        id.data[i] = new T[n];
        for(size_t j = 0; j < n; j++) id.data[i][j] = (T)(i == j);
    }

    return id;
}

template<typename T>
T* Matrix<T>::operator[](size_t index){
    if(index > this -> rows) throw "Index out of bounds";
    return this -> data[index];
}

template<typename T>
size_t Matrix<T>::getRows(){return this -> rows;}

template<typename T>
size_t Matrix<T>::getColumns(){return this -> columns;}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& m){
    if(this -> rows != m.rows || this -> columns != m.columns) throw "Invalid operation: can't sum matrices of different dimensions";

    Matrix result(this -> rows, this -> columns);
    for(size_t i = 0; i < result.rows; i++){
        for(size_t j = 0; j < result.columns; j++){
            result.data[i][j] =this -> data[i][j] + m.data[i][j];
        }
    }

    return result;
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& m){
    if(this -> rows != m.rows || this -> columns != m.columns) throw "Invalid operation: can't sum matrices of different dimensions";

    for(size_t i = 0; i < this -> rows; i++){
        for(size_t j = 0; j < this -> columns; j++){
            (this -> data[i][j]) += m.data[i][j];
        }
    }

    return *(this);
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& m){
    if(this -> rows != m.rows || this -> columns != m.columns) throw "Invalid operation: can't sum matrices of different dimensions";

    Matrix result(this -> rows, this -> columns);
    for(size_t i = 0; i < result.rows; i++){
        for(size_t j = 0; j < result.columns; j++){
            result.data[i][j] = this -> data[i][j] - m.data[i][j];
        }
    }

    return result;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& m){
    if(this -> rows != m.rows || this -> columns != m.columns) throw "Invalid operation: can't sum matrices of different dimensions";

    for(size_t i = 0; i < this -> rows; i++){
        for(size_t j = 0; j < this -> columns; j++){
            this -> data[i][j] -= m.data[i][j];
        }
    }

    return *(this);
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& m){
    if(this -> columns != m.rows) throw "Invalid operation: operation between matrices (A*B) must follow the rule of dimensions A: m*n, B:n*p";
    
    Matrix result = Matrix(this -> rows, m.columns);

    for(size_t i = 0; i < this -> rows; i++){
        for(size_t j = 0; j < m.columns; j++){
            for(size_t k = 0; k < this -> columns; k++){
                result.data[i][j] += this -> data[i][k] * m.data[k][j];
            }
        }
    }

    return result;
}

template<typename T>
Matrix<T> operator*(const T& val, const Matrix<T>& m){
    Matrix<T> result(m);

    for(size_t i = 0; i < result.rows; i++){
        for(size_t j = 0; j < result.columns; j++){
            result.data[i][j] *= val;
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T& val){

    Matrix<T> result(*(this));

    for(size_t i = 0; i < result.rows; i++){
        for(size_t j = 0; j < result.columns; j++){
            result.data[i][j] *= val;
        }
    }

    return result;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const T& val){
    for(size_t i = 0; i < this -> rows; i++){
        for(size_t j = 0; j < this -> columns; j++){
            this -> data[i][j] *= val;
        }
    }

    return *(this);
}

template<typename T>
Matrix<T> Matrix<T>::transpose(const Matrix<T>& m){
    Matrix transposed = Matrix(m.columns, m.rows);

    for(size_t i = 0; i < m.rows; i++){
        for(size_t j = 0; j < m.columns; j++){
            transposed.data[i][j] = m.data[j][i];
        }
    }

    return transposed;
}

template<typename T>
T Matrix<T>::det(const Matrix<T>& m){
    T cumulative_sum = 0;

    if(m.rows != m.columns) throw "Invalid operation: can't compute the determinant of a non square matrix";
    else if(m.rows == 0);
    else if(m.rows == 1) cumulative_sum = m.data[0][0];
    else if(m.rows == 2) cumulative_sum = (m.data[0][0] * m.data[1][1]) - (m.data[0][1] * m.data[1][0]);
    else{
        // implemetation of Laplace expansion (doesnt scale well with bigger matrices: O(n!))

        for(size_t i = 0; i < m.rows; i++){
            if(m.data[i][0] != 0){
                Matrix cofactor;
                cofactor.rows = m.rows - 1;
                cofactor.columns = m.columns - 1;
                cofactor.data = new T*[m.rows - 1];

                bool skip_row = false;
                for(size_t j = 0; j < m.rows; j++){
                    if(i == j) skip_row = true;
                    else cofactor.data[j- skip_row] = &(m.data[j][1]);
                }
                // std::cout << "Cofactor(" << i << "):\n"<< cofactor << "\n";
                cumulative_sum += ((i % 2) ? -1 : 1) * Matrix<T>::det(cofactor) * m.data[i][0];
                // std::cout << "\nPartial sum (" << i << "): " << cumulative_sum << "\n";

                cofactor.data = nullptr;
            }
        }
    }

    // std::cout << "Determinant = "<< cumulative_sum;
    return cumulative_sum;
}

template <typename U>
std::ostream& operator <<(std::ostream &out, const Matrix<U> &m) {
    for(size_t i = 0; i < m.rows; i++){
        for(size_t j = 0; j < m.columns; j++) out << (m.data[i][j]) << "\t";
        out << "\n";
    }

    return out;
}

template <typename T>
Matrix<T>::~Matrix(){
    if(this -> data != nullptr){
        for(size_t i = 0; i < this -> rows; i++) {
            delete[] this -> data[i];
            this -> data[i] = nullptr;
        }
        delete[] this -> data;
        this -> data = nullptr;
    }
}

// ------- Linear algebra functions implementation -------
template<typename T>
void Matrix<T>::swapRows(size_t row1, size_t row2){

    if(row1 >= this -> rows || row2 >= this-> rows) throw "Index out of bounds";

    T* temp = this -> data[row1];
    this -> data[row1] = this -> data[row2];
    this -> data[row2] = temp;
}

template<typename T>
void Matrix<T>::multiplyRow(size_t row, T val){
    if(row >= this -> rows) throw "Index out of bounds";

    size_t m = this -> columns;
    for(int j = 0; j < m; j++) this -> data[row][j] *= val;
}

template<typename T>
void Matrix<T>::addRow(size_t row1, size_t row2, T val){
    if(row1 >= this -> rows || row2 >= this -> rows) throw "Index out of bounds";

    size_t m = this -> columns;
    for(int j = 0; j < m; j++) this -> data[row1][j] += val * this -> data[row2][j];
}

template<typename T>
Matrix<T>* LPU_decomposition(const Matrix<T> &m){

    Matrix<T>* LPU = new Matrix<T>[3];
    LPU[2] = m;

    size_t n = LPU[2].getRows();
    if(n != LPU[2].getColumns()) throw "Invalid format: need a square matrix to perform operation";

    Matrix<T> Id = Matrix<T>::identity(n);
    LPU[0] = Id;
    LPU[1] = Id;

    for(size_t i = 0; i < n - 1; i++){
        if(LPU[2][i][i] == 0){
            size_t j = i + 1;
            bool found_pivot = false;
            while(!found_pivot && j < n){
                if(LPU[2][i][j] != 0){
                    LPU[2].swapRows(i, j);
                    LPU[1].swapRows(i, j);

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
            LPU[0][j][i] = LPU[2][j][i] / LPU[2][i][i];
            LPU[2].addRow(j, i, -LPU[0][j][i]);
        }
    }

    return LPU;
}

template<typename T>
Matrix<T>* QR_decomposition(const Matrix<T> &m){
    Matrix<T>* QR = new Matrix<T>[2];
    QR[0] = Matrix<T>::transpose(m);

    size_t n = QR[0].getRows();
    if(n != QR[0].getColumns()) throw "Invalid format: need a square matrix to perform operation";   

    Vec<T> *u = new Vec<T>[n]; 
    for(size_t i = 0; i < n; i++){
        u[i] = Vec<T>(QR[0][i], n);
        
        for(size_t j = 0; j < i; j++)
            u[i] -= Vec<T>::dot(u[i], u[j]) * u[j];

        u[i] *= sqrt(1 / (double) u[i].getNorm());

        for(size_t j = 0; j < n; j++) QR[0][i][j] = u[i][j];
    }
    QR[1] = QR[0] * m;
    QR[0] = Matrix<T>::transpose(QR[0]);

    return QR; 
}