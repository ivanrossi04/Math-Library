// !V1
// ------- Vector methods implementation -------

template <typename T>
Vec<T>::Vec(){
    this -> val = nullptr;
    this -> dim = 0;
}

template <typename T>
Vec<T>::Vec(const size_t n){
    this -> val = new T[n];
    this -> dim = n;
}

template <typename T>
Vec<T>::Vec(const T arr[],const size_t size){
    if(size) {
        this -> dim = size;
        this -> val = new T[size];
    
        for(size_t i = 0; i < size; i++) this -> val[i] = arr[i];
    } else Vec<T>();
}

template <typename T>
Vec<T>::Vec(std::initializer_list<T> const& arr) {
    this -> dim = arr.size();
    this -> val = new T[this -> dim];
    std::copy(arr.begin(), arr.end(), this -> val);
}

template <typename T>
Vec<T>::Vec(Vec<T> const& v){
    if(v.dim) {
        this -> dim = v.dim;
        this -> val = new T[v.dim];
        for(size_t i = 0; i < v.dim; i++) this -> val[i] = v.val[i];
    } else Vec<T>();
}


// Gru, the Vector Destroyer
template <typename T>
Vec<T>::~Vec(){
    if(this -> val != nullptr) {
        delete[] this -> val;
        this -> val = nullptr;
    }
}

template <typename T>
void Vec<T>::operator=(Vec<T> const& v){
    if(this -> val != nullptr) {
        delete[] this -> val;
        this -> val = nullptr;
    }
    
    if(v.dim) {
        this -> dim = v.dim;
        this -> val = new T[v.dim];
        for(size_t i = 0; i < v.dim; i++) this -> val[i] = v.val[i];
    } else Vec<T>();
}

template <typename T>
void Vec<T>::operator=(std::initializer_list<T> const& arr){
    operator=(Vec<T>(arr));
}

template <typename T>
void Vec<T>::operator=(Vec<T>&& v){
    this -> dim = v.dim;

    if(this -> val != nullptr) delete[] this -> val;
    this -> val = v.val;
    v.val = nullptr;
}

template<typename T>
size_t Vec<T>::getDim(){
    return this -> dim;
}

template<typename T>
T Vec<T>::getNorm(){
    return dot(*(this), *(this));
}

template<typename T>
T& Vec<T>::operator[](size_t index){
    if(index >= this -> dim) throw "Index out of bounds";
    return this -> val[index];
}

template<typename T> bool Vec<T>::operator==(const Vec<T>& rhs){
    if(this -> dim != rhs.dim) return false;
    for(size_t i = 0; i < this -> dim; i++) if(this -> val[i] != rhs.val[i]) return false;
    return true;
}

template<typename T> bool Vec<T>::operator!=(const Vec<T>& rhs){
    return !operator==(rhs);
}

template<typename T> bool Vec<T>::operator<(const Vec<T>& rhs){
    return this -> getNorm() < rhs.getNorm();
}

template<typename T> bool Vec<T>::operator>(const Vec<T>& rhs){
    return this -> getNorm() > rhs.getNorm();
}

template<typename T> bool Vec<T>::operator<=(const Vec<T>& rhs){
    return !operator>(rhs);
}

template<typename T> bool Vec<T>::operator>=(const Vec<T>& rhs){
    return !operator<(rhs);
}

template <typename T>
Vec<T> Vec<T>::operator+(const Vec<T>& vec){

    if((this -> dim) != vec.dim)
        throw "Invalid operation: can't sum vectors of different dimensions";

    Vec<T> result(*this);
    for(size_t i = 0; i < this -> dim; i++)
        result.val[i] += vec.val[i];
    
    return result;
}

template <typename T>
Vec<T>& Vec<T>::operator+=(const Vec<T>& vec){
    if((this -> dim) != vec.dim)
        throw "Invalid operation: can't sum vectors of different dimensions";

    for(size_t i = 0; i < this -> dim; i++)
        this -> val[i] += vec.val[i];
    
    return *this;
}

template <typename T>
Vec<T> Vec<T>::operator-(const Vec<T>& vec){
    if(this -> dim != vec.dim) throw "Invalid operation: can't sum vectors of different dimensions";

    Vec<T> result(*this);
    for(size_t i = 0; i < this -> dim; i++)
        result.val[i] -= vec.val[i];
    
    return result;
}

template <typename T>
Vec<T>& Vec<T>::operator-=(const Vec<T>& vec){
    if(this -> dim != vec.dim) throw "Invalid operation: can't sum vectors of different dimensions";

    for(size_t i = 0; i < this -> dim; i++)
        this -> val[i] -= vec.val[i];
    
    return *this;
}

template <typename T>
Vec<T> Vec<T>::operator*(const T& val){
    Vec<T> result(*this);

    for(size_t i = 0; i < this -> dim; i++)
        result.val[i] *= val;

    return result;
}

template<typename U> Vec<U> operator*(const U& val, const Vec<U>& vec){
    Vec<U> result(vec);

    for(size_t i = 0; i < vec.dim; i++)
        result.val[i] *= val;

    return result;
};

template <typename T>
Vec<T>& Vec<T>::operator*=(const T& val){
    for(size_t i = 0; i < this -> dim; i++)
        this -> val[i] *= val;

    return *this;
}

template <typename T>
T Vec<T>::dot(const Vec<T>& lhs,const Vec<T>& rhs){
    if(lhs.dim != lhs.dim) throw "Invalid operation: can't compute the dot product vectors of different dimensions";
    else if(lhs.dim == 0) throw "Invalid operation: can't compute the dot product vectors of empty vectors";

    T result = lhs.val[0] * rhs.val[0];
    for(size_t i = 1; i < lhs.dim; i++){
        result += lhs.val[i] * rhs.val[i];
    }

    return result;
}

template <typename T>
Vec<T> Vec<T>::cross(const Vec<T>& lhs,const Vec<T>& rhs){
    if(lhs.dim != lhs.dim) throw "Invalid operation: can't compute the cross product vectors of different dimensions";
    else if(lhs.dim != 3) throw "Invalid operation: can't compute the cross product vectors of dim != 3";

    Vec<T> result = Vec<T>({
        lhs.val[1] * rhs.val[2] - lhs.val[2] * rhs.val[1],
        lhs.val[2] * rhs.val[0] - lhs.val[0] * rhs.val[2],
        lhs.val[0] * rhs.val[1] - lhs.val[1] * rhs.val[0]
    });

    return result;
}

template <typename U>
std::ostream& operator <<(std::ostream &out, const Vec<U> &v) {
    for(size_t i = 0; i < v.dim; i++) out << v.val[i] << "\t";

    return out;
}