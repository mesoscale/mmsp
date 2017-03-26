### MMSP.utility.hpp-
void check_boundary(int& x, int x0, int x1, int b0, int b1)
    * changes coordinate x to different values depending on limiting coordinates and boundary conditions
template <int dim, int index, typename T> class target
    * utility class that allows grids to be made with subscripting
    * may be constructed with pointers to data, s0, sx, x0, x1, b0, and b1, or with 
    * operator [] calls check_boundary and returns a target object of <dim-1,index+1,T> with data increased by (x-s0[index+1])*sx[index+1]
    * if dim is 1, operator [] does the same thing but only returns a pointer to the increased data
    * if dim is 0:
        * operator = sets *data to the given value and returns *data
        * operator () returns *data
unsigned long buffer_size(const T& value)
    * returns sizeof(value)
unsigned long to_buffer(const T& value, char* buffer)
    * sends value to the buffer and returns sizeof(value)
unsigned long from_buffer(T& value, const char* buffer)
    * moves a value from the buffer to value and returns sizeof(value)
void read(T& value, std::ifstream& file)
    * reads value from given file
void write(T& value, std::ofstream& file)
    * writes value to given file
int length(const T& value)
    * returns 1
void resize(T& value, int n)
    * does nothing
void copy(T& value, const T& copy)
    * sets value equal to copy
void swap(T& value, T& swap)
    * switches values of value and swap
std::string name(const T& value)
    * returns the type given as a std::string
max(const T& x, const T& y)
    * returns the larger of the two inputs
min(const T& x, const T& y)
    * returns the smaller of the two inputs
T global(T& value, const char* operation)
    * initializes a global value for use in MPI and operates on it depending on the second input, then returns the global value
void print_progress(const int step, const int steps)
    * prints a progress bar and timestamps depending on how far along a process is

### MMSP.vector.hpp
template <typename T> class vector
    * vector class designed to work with MMSP
    * if constructed empty, size is 0 and data is NULL
    * if constructued with another vector, made as a deep copy of the other vector
    * if constructued with an int N, size is set to N and data is set to a dynamic T array of length size
    * if constructed with an int N and a const U& value, size is set to N, data is set to a dynamic array of length size, and each element in data is set to value static_cast to T
    * destructor deletes data and sets to NULL if data is not already NULL
    * operator [] with a given int i returns data[i]
    * operator = sets each value in data to the a given value static_cast to type T, and returns *this
    * operator = with a vector object input deletes data, sets it to NULL, makes the first object a deep copy of the second, and returns *this
    int buffer_size()
        * returns sizeof(size)+size*sizeof(T) (the size of the buffer)
    int to_buffer(char* buffer) const
        * moves size and data to buffer, then returns the buffer size
    int from_buffer(const char* buffer)
        * sets size and data to what is in the buffer, then returns the buffer size
    void write(std::ofstream& file) const
        * writes data to given file
    void read(std::ifstream& file)
        * reads data from given file
    int length() const
        * returns size
    void resize(int N)
        * resizes the vector, leaving last elements empty if the new size is longer or losing last elements if the new size is shorter
    void copy(const vector& v)
        * makes this vector a copy of the given vector
    void swap(vector& v)
        * swaps data and size for two vectors
    void append(const U& value)
        * increases the size of the vector by 1 and places the given value as the last element, static_cast to T
    void append(const vector<U>& v)
        * appends the given vector to this vector, increasing size by the length of the given vector
int buffer_size(const vector<T>& v)
    * returns v.buffer_size()
int to_buffer(const vector<T>& v, char* buffer)
    * returns v.to_buffer(buffer)
int from_buffer(vector<T>& v, const char* buffer)
    * returns v.from_buffer(buffer)
void write(const vector<T>& v, const char* buffer)
    * returns v.write(file)
void read(vector<T>& v, std::ifstream& file)
    * returns v.read(file)
int length(const vector<T>& v)
    * returns v.length
void resize(vector<T>& v, int n)
    * calls v.resize(n)
void copy(vector<T>& v, const vector<T>& w)
    * calls v.copy(w)
void swap(vector<T>& v, vector<T>& w)
    * calls v.swap(w)
void append(vector<T>& v, const U& value)
    * calls v.append(value)
void append(vector<T>& v, const vector<U>& w)
    * calls v.append(w)
std::string name(const vector<T> s)
    * returns "vector:"+name(T())
vector<T> min(const vector<T>& x, const vector<T>& y)
    * returns a vector with each element being the smaller of the corresponding elements in x and y
vector<T> max(const vector<T>& x, const vector<T>& y)
    * returns a vector with each element being the larger of the corresponding elements in x and y
vector<T>& operator+=(vector<T>& x, const vector<U>& y)
    *increases each element in x by its corresponding element in y, then returns x
vector<T> operator+(const vector<T>& x, const vector<U>& y)
    * returns a temporary vector with each element being equal to the corresponding elements in x and y added together
vector<T>& operator-=(vector<T>& x, const vector<U>& y)
    *decreases each element in x by its corresponding element in y, then returns x
vector<T>& operator-(const vector<T>& x, const vector<U>& y)
    * returns a temporary vector with each element being equal to the corresponding elements in y being subtracted from x
vector<T>& operator*=(vector<T>& x, const U& value)
    * multiplies each element in x by the given value
vector<T>& operator*(const U& value, const vector<T>& x)
    * returns a temporary vector consisting of each element in x multiplied by the given value
vector<T>& operator*(const vector<T>&x, const vector<T>& y)
    * returns the dot product of the given vectors
template <int ind, typename T> class target<0,ind,vector<T>>
    * dim=0 specialization for vector class
    * constructor takes pointers to data, s0, sx, x0, x1, b0, and b1 and sets those member variables equal to the pointers
    operator vector<T>&()
        * returns a pointer to data
    operator const vector<T>&() const
        * returns a pointer to data
    T& operator[](int i)
        * returns data->operator[](i) 
    const T& operator[](int i) const
        * returns data->operator[](i) 
    vector<T>& operator=(const T& value) const
        * returns data->operator=(value) 
    vector<T>& operator=(const vector<T>& v) const
        * returns data->operator=(v) 
    vector<T>& operator=(const U& value) const
        * returns data->operator=(value) 
    vector<T>& operator=(const vector<U>& v) const
        * returns data->operator=(v) 
    int buffer_size() const
        * returns data->buffer_size() 
    int to_buffer(char* buffer) const
        * returns data->to_buffer(buffer) 
    int from_buffer(const char* buffer) const
        * returns data->from_buffer(buffer) 
    void write(std::ofstream& file) const
        * calls data->write(file) 
    void read(std::ifstream& file) const
        * calls data->read(file) 
    int length() const
        * returns data->length() 
    int resize(int n) const
        * returns data->resize(n) 
    void copy(const target& t) const
        * calls data->copy(t->data) 
    void swap(const target& t) const
        * calls data->swap(t->data) 
    void append(const U& value) const
        * calls data->append(value) 
    void append(const vector<U>& v) const
        * calls data->append(v)
int buffer_size(const target<0,ind,vector<T>>& v)
    * returns v.buffer_size() 
int to_buffer(const target<0,ind,vector<T>>& v, char* buffer)
    * returns v.to_buffer(buffer) 
int from_bufer(const target<0,ind,vector<T>>& v, const char* buffer)
    * returns v.from_buffer(buffer) 
void write(const target<0,ind,vector<T>>& v,std::ofstream& file)
    * returns v.write(file) 
void read(const target<0,ind,vector<T>>& v,std::ifstream& file)
    * returns v.read(file)
int length(const target<0,ind,vector<T>>& v)
    * returns v.length() 
void resize(const target<0,ind,vector<T>>& v, int n)
    * calls v.resize(n) 
void copy(const target<0,ind,vector<T>>& v, const target<0,ind,vector<T>>& w)
    * calls v.copy(w) 
void swap(const target<0,ind,vector<T>>& v, const target<0,ind,vector<T>>& w)
    * calls v.swap(w)
void append(const target<0,ind,vector<T>>& v, const U& value)
    * calls v.append(value) 
void append(const target<0,ind,vector<T>>& v, const vector<U>& w)
    * calls v.append(w)
std::string name(const target<0,ind,vector<T>>& s)
    * returns "vector:"+name(T())
vector<T> min(const target<0,ind,vector<T>>& x, const target<0,ind,vector<T>>& y)
    * returns min(*(x.data),*(y.data))
vector<T> max(const target<0,ind,vector<T>>& x, const target<0,ind,vector<T>>& y)
    * returns max(*(x.data),*(y.data))
vector<T> operator+=(const target<0,ind,vector<T>>& x, const target<0,ind,vector<U>>& y)
    * returns operator+=(*(x.data),*(y.data))
vector<T> operator+(const target<0,ind,vector<T>>& x, const target<0,ind,vector<U>>& y)
    * returns operator+(*(x.data),*(y.data))
vector<T> operator-=(const target<0,ind,vector<T>>& x, const target<0,ind,vector<U>>& y)
    * returns operator-=(*(x.data),*(y.data))
vector<T> operator-(const target<0,ind,vector<T>>& x, const target<0,ind,vector<U>>& y)
    * returns operator-(*(x.data),*(y.data))
vector<T>& operator*=(target<0,ind,vector<T>>& x, const U& value)
    * returns operator*=(*(x.data),value)
vector<T>& operator*(const U& value, const target<0,ind,vector<T>>& x)
    * returns operator*(value,*(x.data))
bool operator==(const vector<T>& a, const vector<T>& b)
    * returns false if a and b do not have the same length or have nonidentical corresponding elements, or true otherwise

### MMSP.scalar.hpp
template <typename T> class scalar
    * empty constructor does nothing
    * constructor with const T& value sets data to value
    * constructor with const scalar& s sets data to s's data
    * constructor with const U& value or const scalar<U>& s static_casts the given parameter to T and sets data to it
    operator T&()
        * returns data
    operator const T&() const
        * returns data
    scalar& operator=(const T& value)
        * sets data to value and returns *this
    scalar& operator=(const scalar& s)
        * sets data to s's data
    scalar& operator=(const U& value)
        * static_casts value to T, sets data to it, and returns *this
    scalar& operator=(const scalar<U>& s)
        * static_casts s to T, sets data to it, an returns *this
    int buffer_size() const
        * returns the size of T
    int to_buffer(char* buffer) const
        * sends data to buffer, then returns the size of T
    int from_buffer(const char* buffer)
        * sends information from buffer to data, then returns the size of T
    void write(std::ofstream& file) const
        * writes data to file
    void read(std::ifstream& file)
        * reads data from file
    int length() const
        * returns 1
    void resize(int n)
        * does nothing
    void copy(const scalar& s)
        * sets data to s's data
    void swap(scalar& s)
        * switches data from this object with data from s
int buffer_size(const scalar<T>& s)
    * returns s.buffer_size() 
int to_buffer(const scalar<T>& s, char* buffer)
    * returns s.to_buffer(buffer) 
int from_buffer(scalar<T>& s, const char* buffer)
    * returns s.from_buffer(buffer)
void write(const scalar<T>& s, std::ofstream& file)
    * returns s.write(file) 
void read(scalar<T>& s, std::ifstream& file)
    * returns s.read(file) 
int length(const scalar<T>& s)
    * returns s.length()
void resize(scalar<T>& s, int n)
    * calls s.resize(n) 
void copy(scalar<T>& s, const scalar<T>& t)
    * calls s.copy(t) 
void swap(scalar<T>& s, scalar<T>& t)
    * calls s.swap(t) 
std::string name(const scalar<T>& s)
    * returns "scalar:" + name(T())
template <int ind, typename T> class target<0,ind,scalar<T>>
    * dim=0 specialization for scalar class
    * constructor takes pointers to scalar<T> data, const ints s0, sx, x0, x1, b0, b1, and sets those variables in the objects to the pointers
    operator T&()
        * returns pointer to data
    operator const T&() const
        * returns pointer to data
    scalar<T>& operator=(const T& value) const
        * returns data->operator=(value) 
    scalar<T>& operator=(const scalar<T>& s)
        * returns data->operator=(s) 
    scalar<T>& operator=(const U& value) const
        * returns data->operator=(value) 
    scalar<T>& operator=(const scalar<T>& s) const
        * returns data->operator=(s) 
    int buffer_size() const
        * returns data->buffer_size() 
    int to_buffer(char* buffer) const
        * returns data->to_buffer(buffer) 
    int from_buffer(const char* buffer) const
        * returns data->from_buffer(buffer) 
    void write(std::ofstream& file) const
        * calls data->write(file) 
    void read(std::ifstream& file) const
        * calls data->read(file) 
    int length() const
        * returns data->length() 
    int resize(int n) const
        * returns data->resize(n) 
    void copy(const target& t) const
        * calls data->copy(t->data) 
    void swap(const starget& t) const
        * calls data->swap(t->data) 
int buffer_size(const target<0,ind,scalar<T>>& s)
    * returns s.buffer_size() 
int to_buffer(const target<0,ind,scalar<T>>& s, char* buffer)
    * returns s.to_buffer(buffer) 
int from_buffer(const target<0,ind,scalar<T>>& s, const char* buffer)
    * returns s.from_buffer(buffer) 
void write(const target<0,ind,scalar<T>>& s, std::ofstream& file)
    * returns s.write(file) 
void read(const target<0,ind,scalar<T>>& s, std::ifstream& file)
    * returns s.read(file)
int length(const target<0,ind,scalar<T>>& s)
    * returns s.length() 
void resize(const target<0,ind,scalar<T>>& s, int n)
    * calls s.resize(n)
void copy(const target<0,ind,scalar<T>>& s, const target<0,ind,scalar<T>>& t)
    * calls s.copy(t) 
void swap(const target<0,ind,scalar<T>>& s, const target<0,ind,scalar<T>>& t)
    * calls s.swap(t)
std::string name(const target<0,ind,scalar<T>>& s)
    * returns "scalar:" + name(T())

### MMSP.sparse.hpp-
template <typename T> struct item
    * holds int index and T value
template <typename T> class sparse:
    * class for sparsely populated arrays, done by an array of item objects
    * empty constructor sets size to 0 and data to NULL
    * constructor with a sparse object input sets size and data to that of the input sparse object
    * destructor deletes data and sets to NULL if it is not already
    sparse& operator=(const sparse& x)
        * sets size and data to that of x, and returns *this
    sparse& operator=(const sparse<U>& x)
        * sets size and data to that of x, static casting each value in data, then returns *this
    T& set(int index)
        * creates a new item in the object with the given index and a value of 0
    T operator[](int index) const
        * returns the value associated with the index in the object, or 0 if there is none
    unsigned int grain_id() const
        * returns the index associated with the highest value in the sparse object
    double getMagPhi() const
        * returns the square root of the sum of the squares of each value in the object
    int index(int i) const
        * returns the index associated with the i-th item in the object
    T value(int i) const
        * returns the value associated with the i-th item in the object
    int buffer_size() const
        * returns sizeof(size) + size * sizeof(item<T>) (the size of the buffer)
    int to_buffer(char* buffer)
        * moves data into buffer, then returns the size of the buffer
    int from_buffer(char* buffer)
        * moves data from the buffer, then returns the size of the buffer
    void write(std::ofstream& file)
        * writes size and data to file
    void read(std::ifstream file)
        * reads size and data from file
    int length() const
        * returns size
    void resize(int n)
        * does nothing
    void copy(const sparse& s)
        * makes this object a deep copy of s
    void swap(sparse& s)
        * switches size and data of this object and s
int buffer_size(const sparse<T>& s)
    * returns s.buffer_size()
int to_buffer(const sparse<T>& s, char* buffer)
    * returns s.to_buffer(buffer)
int from_buffer(sparse<T>& s, const char* buffer)
    * returns s.from_buffer(buffer)
void write(const sparse<T>& s, std::ofstream& file)
    * returns s.write(file)
void read(const sparse<T>& s, std::ifstream& file)
    * returns s.read(file)
int length(const sparse<T>& s)
    * returns s.length()
void resize(sparse<T>& s, int n)
    * calls s.resize(n)
void copy(sparse<T>& s, const sparse<T>& t)
    * calls s.copy(t)
void swap(sparse<T>& s, sparse<T>& t)
    * calls s.swap(t)
T& set(sparse<T>& s, int index)
    * calls s.set(index)
int index(const sparse<T>& s, int i)
    * returns s.index(i)
T value(const sparse<T>& s, int i)
    * returns s.value(i)
std::string name(const sparse<T>& s)
    * returns "sparse:" + name(T())
sparse<T> max(const sparse<T>& x, const sparse<T>& y)
    * returns a sparse object filled with the larger values from x or y at the proper indices
sparse<T> min(const sparse<T>& x, const sparse<T>& y)
    * returns a sparse object filled with the smaller values from x or y at the proper indices    
sparse<T>& operator+=(sparse<T>&x, const sparse<U>& y)
    * increases each item in x by the value of its corresponding item in y, static_casting to do so
sparse<T>& operator-=(sparse<T>& x, const sparse<U>& y)
    * decreases each item in x by the value of its corresponding item in y, static_casting to do so
sparse<T> operator-(const sparse<T>& x, const sparse<U>& y)
    * returns a sparse object populated with the items in y subtracted from the items in x
sparse<T>& operator*=(sparse<T>& x, const U& value)
    * permanently multiplies each item in x by value
sparse<T> operator*(const U& value, const sparse<T>& x)
    * returns a sparse object equal to each item in x multiplied by value
template <int ind, typename T> class target<0,ind,sparse<T>>
    * specialization for target class dim=0
    * constructor takes pointers to data, s0, sx, x0, x1, b0, b1, and sets member variables to these
    operator sparse<T>&()
        * returns *data
    operator const sparse<T>&() const
        * returns *data
    T operator[](int i)
        * returns data->operator[](i)
    sparse<T>& operator=(const sparse<T>& s) const
        * returns data->operator=(s)
    sparse<T>& operator=(const sparse<U>& s) const
        * returns data->operator=(s)
    int buffer_size() const
        * returns data->buffer_size()
    int to_buffer(char* buffer) const
        * returns data->to_buffer(buffer)
    int from_buffer(const char* buffer) const
        * returns data->from_buffer(buffer)
    void write(std::ofstream& file) const
        * calls data->write(file)
    void read(std::ofstream& file) const
        * calls data->read(file)
    int length() const
        * returns data->length
    int resize(int n) const
        * returns data->resize(n)
    void copy(const target& t) const
        * returns data->copy(t->data)
    void swap(const target& t) const
        * returns data->swap(t->data)
    T& set(int index) const
        * returns data->set(index)
    int index(int i) const
        * returns data->index(i)
    T value(int i) const
        * returns data->value(i)
T& set(const target<0,ind,sparse<T>>& s, int index)
    * returns s.set(index)
int index(const target<0,ind,sparse<T>>& s,int i)
    * returns s.index(i)
T value(const target<0,ind,sparse<T>>& s,int i)
    * returns s.value(i)
int buffer_size(const target<0,ind,sparse<T>>& s, int index)
    * returns s.buffer_size()
int to_buffer(const target<0,ind,sparse<T>>& s, char* buffer)
    * returns s.to_buffer(buffer)
int from_buffer(const target<0,ind,sparse<T>>& s, const char* buffer)
    * returns s.from_buffer(buffer)
void write(const target<0,ind,sparse<T>>& s, std::ofstream& file)
    * returns s.write(file)
void read(const target<0,ind,sparse<T>>& s, std::ifstream& file)
    * returns s.read(file)
int length(const target<0,ind,sparse<T>>& s)
    * returns s.length()
void resize(const target<0,ind,sparse<T>>& s, int n)
    * calls s.resize(n)
void copy(const target<0,ind,sparse<T>>& s, const target<0,ind,sparse<T>>& t)
    * calls s.copy(t)
void swap(const target<0,ind,sparse<T>>& s, const target<0,ind,sparse<T>>& t)
    * calls s.swap(t)
std::string name(const target<0,ind,sparse<T>>& s)
    * returns "sparse:" + name(T())
sparse<T>& min(target<0,ind,sparse<T>> c, const target<0,ind,sparse<T>>& y)
    * returns min(*(x.data),*(y.data))
sparse<T>& max(target<0,ind,sparse<T>> c, const target<0,ind,sparse<T>>& y)
    * returns max(*(x.data),*(y.data))
sparse<T>& operator+=(target<0,ind,sparse<T>> c, const target<0,ind,sparse<U>>& y)
    * returns operator+=(*(x.data),*(y.data))
sparse<T> operator+(target<0,ind,sparse<T>> c, const target<0,ind,sparse<U>>& y)
    * returns operator+(*(x.data),*(y.data))
sparse<T>& operator-=(target<0,ind,sparse<T>> c, const target<0,ind,sparse<U>>& y)
    * returns operator-=(*(x.data),*(y.data))
sparse<T> operator-(target<0,ind,sparse<T>> c, const target<0,ind,sparse<U>>& y)
    * returns operator-(*(x.data),*(y.data))
sparse<T>& operator*=(target<0,ind,sparse<T>> c, const U& value)
    * returns operator-=(*(x.data),value)
sparse<T> operator*(const U& value, const target<0,ind,sparse<T>>& x)
    * returns operator*(value,*(x.data))
bool operator==(const sparse<T>& a, const sparse<T>& b)
    * returns true if all indices of both objects contain the same values, and false otherwise
    
### MMSP.grid.hpp
int nodes(const grid<dim,T>& GRID)
    * returns nodes(GRID)
int fields(const grid<dim,T>& GRID)
    * returns fields(GRID)
int ghosts(const grid<dim,T>& GRID)
    * returns ghosts(GRID)
int g0(const grid<dim,T>& GRID, int i)
    * returns g0(GRID, i)
int g1(const grid<dim,T>& GRID, int i)
    * returns g1(GRID, i)
int b0(const grid<dim,T>& GRID, int i)
    * returns b0(GRID, i)
int b1(const grid<dim,T>& GRID, int i)
    * returns b1(GRID, i)
int& b0(const grid<dim,T>& GRID, int i)
    * returns b0(GRID, i)
int& b1(const grid<dim,T>& GRID, int i)
    * returns b1(GRID, i)
int x0(const grid<dim,T>& GRID, int i)
    * returns x0(GRID, i)
int x1(const grid<dim,T>& GRID, int i)
    * returns x1(GRID, i)
int xmin(const grid<dim,T>& GRID, int i)
    * returns xmin(GRID, i)
int xmax(const grid<dim,T>& GRID, int i)
    * returns xmax(GRID, i)
int xlength(const grid<dim,T>& GRID, int i)
    * returns xlength(GRID, i)
double dx(const grid<dim,T>& GRID, int i)
    * returns dx(GRID, i)
double& dx(grid<dim,T>& GRID, int i)
    * returns dx(GRID, i)
int x0(const grid<dim,T>& GRID)
    * returns x0(GRID)
int x1(const grid<dim,T>& GRID)
    * returns x1(GRID)
int xmin(const grid<dim,T>& GRID)
    * returns xmin(GRID)
int xmax(const grid<dim,T>& GRID)
    * returns xmax(GRID)
int xlength(const grid<dim,T>& GRID)
    * returns xlength(GRID)
double dx(const grid<dim,T>& GRID)
    * returns dx(GRID)
double& dx(grid<dim,T>& GRID)
    * returns dx(GRID)
int y0(const grid<dim,T>& GRID)
    * returns y0(GRID)
int y1(const grid<dim,T>& GRID)
    * returns y1(GRID)
int ymin(const grid<dim,T>& GRID)
    * returns ymin(GRID)
int ymax(const grid<dim,T>& GRID)
    * returns ymax(GRID)
int ylength(const grid<dim,T>& GRID)
    * returns ylength(GRID)
double dy(const grid<dim,T>& GRID)
    * returns dy(GRID)
double& dy(grid<dim,T>& GRID)
    * returns dy(GRID)
int z0(const grid<dim,T>& GRID)
    * returns z0(GRID)
int z1(const grid<dim,T>& GRID)
    * returns z1(GRID)
int zmin(const grid<dim,T>& GRID)
    * returns zmin(GRID)
int zmax(const grid<dim,T>& GRID)
    * returns zmax(GRID)
int zlength(const grid<dim,T>& GRID)
    * returns zlength(GRID)
double dz(const grid<dim,T>& GRID)
    * returns dz(GRID)
double& dz(grid<dim,T>& GRID)
    * returns dz(GRID)
template <int dim, typename T> class grid
    * central class of MMSP, allows for calculation of mathematical functions of arrays of 2 or 3 dimensions with nearly identical function calls
    * constructor with int FIELDS, int min[dim], int max[dim], int GHOSTS (with default of 1) and bool SINGLE (with default of false) sets given variables in the grid, sets each index in g0 and g1 to those in min and max, respectively, and calls setup(SINGLE)
    * constructor with int FIELDS first and a list of arguments after will set fields to FIELDS and use each other argument to populate g0 and g1
    * constructor with a const grid object will make this object a deep copy of that object
    * constructor with a const grid object of type u and int FIELDS will set fields to FIELDS and make this grid object a deep copy of the given object
    * constructor with const char* filename and int GHOSTS (with default of 1) will set data to NULL and call input(filename, GHOSTS) to read data from the file
    void setup(bool SINGLE=false)
        *sets grid default values and does the base work to allow parallel computation if the MPI version is being used
    * destructor deletes data and sets to NULL
    grid& operator=(const U& value)
        * sets each entry in data to value static_cast to T
    grid& operator=(const grid<dim,U>& GRID)
        * sets each entry in data to its corresponding entry in GRID.data static_cast to T
    grid& operator+=(const grid<dim,U)& GRID)
        * increases each entry in data by its corresponding entry in GRID.data static_cast to T
    grid& operator-=(const grid<dim,U)& GRID)
        * decreases each entry in data by its correspondng entry in GRID.data static_cast to T
    target<dim-1,0,T> operator [](int x) const
        * calls check_boundary(x,x0[0],x1[0],b0[0],b1[0]) and returns a target<dim-1,0,T> object initialized with (data+(x-s0[0])*sx[0],sx,x0,x1,b0,b1)
    T& operator()(MMSP::vector<int> x) const
        * returns a pointer to data after calling check_boundary(x[i],x0[i],x1[i],b0[i],b1[i]) for each dimension and increasing the pointer by (x[i]-s0[i])*sx[i] for each dimension
    T& operator()(int i) const
        * returns data[i], changing the value of i depending on local coordinate limits if using the MPI version
    MMSP::vector<int> position(int i) const
        * returns a vector stating the position in terms of each dimension for index i in the grid
    void ghostswap()
        * parallelization function, sends computations from one processor to another
    void ghostswap(const int sublattice)
        * parallelization function for Monte Carlo muniation reduction
    unsigned long buffer_size_save(const int min[dim], const int max[dim]) const
        * returns buffer_size_save(data,0,min,max)
    unsigned long buffer_size_save(T* p, int i, const int min[dim], const int max[dim]) const
        * returns the buffer size for the given pointer, accounting for each dimension above i, incrementing between min and max by increments of 2
    unsigned long to_buffer_save(char* buffer, const int min[dim], const int max[dim]) const
        * returns to_buffer_save(buffer, data, 0, min, max)
    unsigned long to_buffer_save(char* buffer, T* p, int i, const int min[dim], const int max[dim]) const
        * does the same as buffer_size_save, except sends a modified pointer to buffer using MMSP::to_buffer
    unsigned long from_buffer_save(char* buffer, const int min[dim], const int max[dim])
        * returns from_buffer_save(buffer, data, 0, min, max)
    unsigned long from_buffer_save(char* buffer, T* p, int i, const int min[dim], const int max[dim]) const
        * does the same as buffer_size_save, except receives a modified pointer to buffer using MMSP::from_buffer
    unsigned long buffer_size() const
        * returns buffer_size(x0,x1)
    unsigned long buffer_size(const int min[dim], const int max[dim]) const
        * returns buffer_size(data, 0, min, max)
    unsigned long buffer_size(T* p, int i, const int min[dim], const int max[dim]) const
        * returns the buffer size for the given pointer, accounting for each dimension above i, incrementing between min and max by increments of 1
    unsigned long to_buffer(char* buffer) const
        * returns to_buffer(buffer, x0, x1)
    unsigned long to_buffer(char* buffer, const int min[dim], const int max[dim]) const
        * returns to_buffer(buffer, data, 0 ,min, max)
    unsigned long to_buffer(char* buffer, T* p, int i, const int min[dim], const int max[dim]) const
        * does the same as buffer_size but also populates the buffer with a modified pointer p using MMSP::to_buffer 
    unsigned long from_buffer(char* buffer) const
        * returns from_buffer(buffer, x0, x1)
    unsigned long from_buffer(char* buffer, const int min[dim], const int max[dim]) const
        * returns from_buffer(buffer, data, 0 ,min, max)
    unsigned long from_buffer(char* buffer, T* p, int i, const int min[dim], const int max[dim]) const
        * does the same as buffer_size but also populates a modified pointer from the buffer using MMSP::to_buffer
    void input(const char* filename, int GHOSTS=1, int SINGLE=false)
        * sets the grid object to an unpopulated grid of the parameters specified in the given file
    void read(std::ifstream& file, int GHOSTS=1)
        * copies data from the given file into the grid where they overlap spatially
    void output(const char* filename) const
        * writes grid parameters and data to file, with versions to work whether or not the MPI version is being used
    unsigned long write_buffer(char*& buf) const
        * writes grid parameters and data to buf
    friend int nodes(const grid& GRID)
        * returns GRID.nodes
    friend int fields(const grid& GRID)
        * returns GRID.fields
    friend int ghosts(const grid& GRID)
        * returns GRID.ghosts
    friend int g0(const grid& GRID, int i)
        * returns GRID.g0[i]
    friend int g1(const grid& GRID, int i)
        * returns GRID.g1[i]
    friend int b0(const grid& GRID, int i)
        * returns GRID.b0[i]
    friend int b1(const grid& GRID, int i)
        * returns GRID.b1[i]
    friend int& b0(grid& GRID, int i)
        * returns GRID.b0[i]
    friend int& b1(grid& GRID, int i)
        * returns GRID.b1[i]
    friend int x0(const grid& GRID, int i)
        * returns GRID.x0[i]
    friend int x1(const grid& GRID, int i)
        * returns GRID.x1[i]
    friend int xmin(const grid& GRID, int i)
        * returns GRID.x0[i]
    friend int xmax(const grid& GRID, int i)
        * returns GRID.x1[i]
    friend int xlength(const grid& GRID, int i)
        * returns GRID.x1[i]-GRID.x0[i]
    friend double dx(const grid& GRID, int i)
        * returns GRID.dx[i]
    friend double& dx(grid& GRID, int i)
        * returns GRID.dx[i]
    friend int N0(const grid& GRID, int i)
        * returns GRID.n0[i]
    friend int N1(const grid& GRID, int i)
        * returns GRID.n1[i]
    friend int P0(const grid& GRID, int i)
        * returns GRID.p0[i]
    friend int P1(const grid& GRID, int i)
        * returns GRID.p1[i]
    friend int sp(const grid& GRID, int i)
        * returns GRID.sp[i]
    friend int x0(const grid& GRID)
        * returns GRID.x0[0]
    friend int x1(const grid& GRID)
        * returns GRID.x1[0]
    friend int xmin(const grid& GRID)
        * returns GRID.x0[0]
    friend int xmax(const grid& GRID)
        * returns GRID.x1[0]
    friend int xlength(const grid& GRID)
        * returns GRID.x1[0]-GRID.x0[0]
    friend double dx(const grid& GRID)
        * returns GRID.dx[0]
    friend double& dx(grid& GRID)
        * returns GRID.dx[0]
    friend int y0(const grid& GRID)
        * returns GRID.x0[1]
    friend int y1(const grid& GRID)
        * returns GRID.x1[1]
    friend int ymin(const grid& GRID)
        * returns GRID.x0[1]
    friend int ymax(const grid& GRID)
        * returns GRID.x1[1]
    friend int ylength(const grid& GRID)
        * returns GRID.x1[1]-GRID.x0[1]
    friend double dy(const grid& GRID)
        * returns GRID.dx[1]
    friend double& dy(grid& GRID)
        * returns GRID.dx[1]
    friend int z0(const grid& GRID)
        * returns GRID.x0[2]
    friend int z1(const grid& GRID)
        * returns GRID.x1[2]
    friend int zmin(const grid& GRID)
        * returns GRID.x0[2]
    friend int zmax(const grid& GRID)
        * returns GRID.x1[2]
    friend int zlength(const grid& GRID)
        * returns GRID.x1[2]-GRID.x0[2]
    friend double dz(const grid& GRID)
        * returns GRID.dx[2]
    friend double& dz(grid& GRID)
        * returns GRID.dx[2]
    void swap(grid& GRID)
        * exchanges data and all parameters between this object and GRID
    void copy(const grid& GRID)
        * makes this object a deep copy of GRID
T laplacian(const grid<dim,T>& GRID, const vector<int>& x)
    * returns the laplacian at the point specified by x for a single field grid
vector<T> laplacian(const grid<dim,vector<T>>& GRID, const vector<int>& x)
    * returns the laplacian for each field in the grid as a vector at the point specified by x
T laplacian(const grid<dim,vector<T>>& GRID, const vector<int>& x, const int field)
    * returns the laplacian for the specified field at the point specified by x
sparse<T> laplacian(const grid<dim,sparse<T>>& GRID, const vector<int>& x)
    * returns the laplacian at the point specified by x as a sparse object
T laplacian(const grid<dim,T>& GRID, int i)
    * returns the laplacian at GRID.position(i) for a single field grid
vector<T> laplacian(const grid<dim,vector<T>>& GRID, int i)
    * returns the laplacian at GRID.position(i) for each field in the grid as a vector 
T laplacian(const grid<dim, vector<T>>& GRID, int i, int f)
    * returns the laplacian at GRID.position(i) for field f
sparse<T> laplacian(const grid<dim, sparse<T>>& GRID, int i)
    * returns the laplacian at GRID.position(i) as a sparse object
vector<T> gradient(const grid<dim, T>& GRID, const vector<int>& x)
    * returns the gradient at the point specified by x
vector<T> grad(const grid<dim,T>& GRID, const vector<int>& x)
    * returns gradient(GRID, x)
T divergence(const grid<dim,T>& GRID, const vector<int>& x)
    * returns the divergence for a single field grid at the point specified by x
vector<T> divergence(const grid<dim,vector<T>>& GRID, const vector<int>& x)
    * returns the divergence for each field in the grid at the point specified by x
T div(const grid<dim,T>& GRID, const vector<int>& x)
    * returns divergence(GRID, x)
vector<T> div(const grid<dim,vector<T>>& GRID, const vector<int>& x)
    * returns divergence(GRID, x)
MMSP::vector<int> position(const grid<dim,T>& GRID, int index)
    * returns GRID.position(index)
void ghostswap(grid<dim,T>& GRID)
    * calls GRID.ghostswap()
void ghostswap(grid<dim,T>& GRID, const int sublattice)
    * calls GRID.ghostswap(sublattice)
unsigned long buffer_size(const grid<dim,T>& GRID)
    * returns GRID.buffer_size()
unsigned long to_buffer(const grid<dim,T>& GRID, char* buffer)
    * returns GRID.to_buffer(buffer)
unsigned long from_buffer(grid<dim,T>& GRID, const char* buffer)
    * returns GRID.from_buffer(buffer)
void read(grid<dim,T>& GRID, std::ifstream& file)
    * calls GRID.read(file)
void write(const grid<dim,T>& GRID, std::ifstream& file)
    * calls GRID.write(file)
void input(grid<dim,T>& GRID, const char* filename, int GHOSTS=1, int SINGLE=false)
    * calls GRID.input(filename, GHOSTS, SINGLE)
void output(const grid<dim,T>& GRID, const char* filename)
    * calls GRID.output(filename)
unsigned long write_buffer(const grid<dim,T>& GRID, char*& buf)
    * calls GRID.write_buffer(buf)
int length(const grid<dim,T>& GRID)
    * returns nodes(GRID)
void resize(int n, grid<dim,T>& GRID)
    * does nothing
void swap(grid<dim,T>& GRID1, grid<dim,T>& GRID2)
    * calls GRID1.swap(GRID2)
void copy(grid<dim,T>& GRID1, grid<dim,T>& GRID2)
    * calls GRID2.copy(GRID1)
std::string name(const grid<dim,T>& GRID)
    * returns "grid:" + name(T())
