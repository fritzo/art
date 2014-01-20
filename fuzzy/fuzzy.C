
#include <map>
#include <fstream>
#include <iostream>
#include <cstdlib>

//fibs = 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040, 1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817, 39088169, 63245986, 102334155
#define FIB01 1U
#define FIB02 1U
#define FIB03 2U
#define FIB04 3U
#define FIB05 5U
#define FIB06 8U
#define FIB07 13U
#define FIB08 21U
#define FIB09 34U
#define FIB10 55U
#define FIB11 89U
#define FIB12 144U
#define FIB13 233U
#define FIB14 377U
#define FIB15 610U
#define FIB16 987U
#define FIB17 1597U
#define FIB18 2584U
#define FIB19 4181U
#define FIB20 6765U
#define FIB21 10946U
#define FIB22 17711U

#define X_MAX 1021
#define Y_MAX 1019
#define X_POS 1
#define Y_POS 2

template<unsigned MODULUS>
class Num
{
    typedef Num<MODULUS> MyType;
    unsigned m_num;
public:
    Num (unsigned n) : m_num(n % MODULUS) {}

    unsigned  to_int() const { return m_num; }
    unsigned& to_int()       { return m_num; }

    bool operator != (const MyType& other) const { return m_num != other.m_num; }
    bool operator < (const MyType& other) const { return m_num < other.m_num; }
    MyType operator + (const MyType& other) const { return m_num + other.m_num; }
    MyType operator - (const MyType& other) const { return m_num - other.m_num + MODULUS; }
    MyType operator * (const MyType& other) const { return m_num * other.m_num; }
};

class Point
{
    typedef Num<X_MAX> Real;
    typedef Num<Y_MAX>  Comp;
    Real m_real;
    Comp m_comp;
public:
    Point () : m_real(0), m_comp(0) {}
    Point (Real x, Comp y) : m_real(x), m_comp(y) {}

    unsigned  real () const { return m_real.to_int(); }
    unsigned  comp () const { return m_comp.to_int(); }
    unsigned& real ()       { return m_real.to_int(); }
    unsigned& comp ()       { return m_comp.to_int(); }

    //sorting operations, for data structures, not math
    bool operator!= (const Point& other) const
    { return (m_real!=other.m_real) or (m_comp!=other.m_comp); }
    bool operator< (const Point& other) const
    { return (m_real!=other.m_real) ? (m_real<other.m_real) : (m_comp<other.m_comp); }

    //commutative operations
    Point operator+ (const Point& other) const;
    Point operator* (const Point& other) const;
};
Point Point::operator+ (const Point& other) const
{
    return Point(m_real + other.m_real,
                 m_comp + other.m_comp);
}
Point Point::operator* (const Point& other) const
{
    Real x1 = m_real.to_int() * other.m_real.to_int();
    Real x2 = m_comp.to_int() * other.m_comp.to_int();
    Comp y1 = m_real.to_int() * other.m_comp.to_int();
    Comp y2 = m_comp.to_int() * other.m_real.to_int();
    return Point(x1 - x2, y1 + y2);
}

class SparsePMF
{
    typedef std::map<Point, float> PMF;
    PMF m_pmf;
public:
    SparsePMF () {}
    SparsePMF (Point x) { m_pmf[x] = 1.0; }

    unsigned size () const { return m_pmf.size(); }

    //output
    void writeToFile() const;

    //measure operations
    float  operator [] (Point x) const;
    float& operator [] (Point x) { return m_pmf[x]; }
    SparsePMF& operator *= (float scale);
    SparsePMF& operator += (const SparsePMF& other);

    //algebraic operations
    friend SparsePMF product (const SparsePMF& lhs, const SparsePMF& rhs);
    friend SparsePMF sum     (const SparsePMF& lhs, const SparsePMF& rhs);

    //iteration
    typedef PMF::iterator iterator;
    iterator begin () { return m_pmf.begin(); }
    iterator end   () { return m_pmf.end(); }

    //const iteration
    typedef PMF::const_iterator const_iterator;
    const_iterator begin () const { return m_pmf.begin(); }
    const_iterator end   () const { return m_pmf.end(); }
};
void SparsePMF::writeToFile() const
{
    std::ofstream f("fuzzy.text");
    for (const_iterator iter = begin(); iter != end(); ++iter) {
        Point x = iter->first;
        float p = iter->second;
        f << x.real() << '\t' << x.comp() << '\t' << p << '\n';
    }
    f.close();
}
float SparsePMF::operator[] (Point x) const
{
    const_iterator iter = m_pmf.find(x);
    if (iter == end()) return 0;
    return iter->second;
}
SparsePMF& SparsePMF::operator*= (float scale)
{
    for (iterator iter = begin(); iter != end(); ++iter) {
        iter->second *= scale;
    }
}
SparsePMF& SparsePMF::operator+= (const SparsePMF& other)
{
    for (const_iterator iter = other.begin(); iter != other.end(); ++iter) {
        m_pmf[iter->first] += iter->second;
    }
}
SparsePMF product (const SparsePMF& lhs, const SparsePMF& rhs)
{//pointwise multiplication
    SparsePMF result;
    SparsePMF::const_iterator iter_lhs, begin_lhs=lhs.begin(), end_lhs=lhs.end();
    SparsePMF::const_iterator iter_rhs, begin_rhs=rhs.begin(), end_rhs=rhs.end();
    for (iter_lhs = begin_lhs; iter_lhs != end_lhs; ++iter_lhs) {
        Point x1 = iter_lhs->first;
        float p1 = iter_lhs->second;
        for (iter_rhs = begin_rhs; iter_rhs != end_rhs; ++iter_rhs) {
            Point x2 = iter_rhs->first;
            float p2 = iter_rhs->second;

            Point x = x1 * x2;
            float p = p1 * p2;

            result[x] += p;
        }
    }
    return result;
}
SparsePMF sum (const SparsePMF& lhs, const SparsePMF& rhs)
{//pointwise addition
    SparsePMF result;
    SparsePMF::const_iterator iter_lhs, begin_lhs=lhs.begin(), end_lhs=lhs.end();
    SparsePMF::const_iterator iter_rhs, begin_rhs=rhs.begin(), end_rhs=rhs.end();
    for (iter_lhs = lhs.begin(); iter_lhs != lhs.end(); ++iter_lhs) {
        Point x1 = iter_lhs->first;
        float p1 = iter_lhs->second;
        for (iter_rhs = rhs.begin(); iter_rhs != rhs.end(); ++iter_rhs) {
            Point x2 = iter_rhs->first;
            float p2 = iter_rhs->second;

            Point x = x1 + x2;
            float p = p1 * p2;

            result[x] += p;
        }
    }
    return result;
}

class DensePMF
{
    float *m_pmf;
    enum { MAX = X_MAX * Y_MAX };
public:
    DensePMF ()         { m_pmf = new float[MAX];
                          for (int i=0; i<MAX; ++i) m_pmf[i] = 0.0f; }
    DensePMF (Point x)  { m_pmf = new float[MAX];
                          for (int i=0; i<MAX; ++i) m_pmf[i] = 0.0f;
                          get(x) = 1.0; }
    ~DensePMF( ) { delete[] m_pmf; }
    DensePMF& operator= (const DensePMF& other);

    unsigned size () const;
    float mass () const;

    //output
    void writeToFile() const;

    //measure operations
    float  get (unsigned i, unsigned j) const { return m_pmf[i + X_MAX*j]; }
    float& get (unsigned i, unsigned j)       { return m_pmf[i + X_MAX*j]; }
    float  get (Point x) const { return get(x.real(), x.comp()); }
    float& get (Point x)       { return get(x.real(), x.comp()); }
    float  operator [] (Point x) const { return get(x); }
    float& operator [] (Point x)       { return get(x); }
    DensePMF& operator *= (float scale);
    DensePMF& operator += (const DensePMF& other);
    void prune (float thresh);

    //algebraic operations
    friend DensePMF product_and_sum (const DensePMF& other);

    //iteration
    class iterator
    {
        DensePMF* pmf;
        unsigned i, j;
    public:
        iterator () {}
        iterator (DensePMF* _pmf, unsigned _i, unsigned _j) : pmf(_pmf), i(_i), j(_j) {}
        bool finished () const { return i==X_MAX; }
        bool operator != (const iterator& other) { return (i!=other.i) or (j!=other.j); }
        iterator& operator ++ () { ++j; if (j == Y_MAX) { j=0; ++i; } return *this; }
        Point  first () { return Point(i,j); }
        float& second() { return pmf->get(i,j); }
    };
    iterator begin () { return iterator(this, 0, 0); }
    iterator end   () { return iterator(this, X_MAX, 0); }

    //const iteration
    class const_iterator
    {
        const DensePMF* pmf;
        unsigned i, j;
    public:
        const_iterator () {}
        const_iterator (const DensePMF* _pmf, unsigned _i, unsigned _j)
            : pmf(_pmf), i(_i), j(_j) {}
        bool finished () const { return i==X_MAX; }
        bool operator != (const const_iterator& other) { return (i!=other.i) or (j!=other.j); }
        const_iterator& operator ++ () { ++j; if (j == Y_MAX) { j=0; ++i; } return *this; }
        Point first () { return Point(i,j); }
        float second() { return pmf->get(i,j); }
    };
    const_iterator begin () const { return const_iterator(this, 0, 0); }
    const_iterator end   () const { return const_iterator(this, X_MAX, 0); }
};
DensePMF& DensePMF::operator= (const DensePMF& other)
{
    for (int i=0; i<MAX; ++i) m_pmf[i] = other.m_pmf[i];
}
unsigned DensePMF::size () const
{
    unsigned result=0;
    for (int i=0; i<MAX; ++i) {
        if (m_pmf[i] > 0) ++result;
    }
    return result;
}
float DensePMF::mass () const
{
    double result=0.0;
    for (int i=0; i<MAX; ++i) result += m_pmf[i];
    return result;
}
void DensePMF::writeToFile() const
{
    std::ofstream f("fuzzy.text");
    for (const_iterator iter = begin(); iter != end(); ++iter) {
        float p = iter.second();
        if (p == 0) continue;
        Point x = iter.first();
        float re = x.real(); re /= X_MAX; if (re > 0.5) re -= 1.0;
        float im = x.comp(); im /= Y_MAX; if (im > 0.5) im -= 1.0;
        f << im << '\t' << re << '\t' << p << '\n';
    }
    f.close();
}
DensePMF& DensePMF::operator*= (float scale)
{
    for (int i=0; i<MAX; ++i) m_pmf[i] *= scale;
}
DensePMF& DensePMF::operator+= (const DensePMF& other)
{
    for (int i=0; i<MAX; ++i) m_pmf[i] += other.m_pmf[i];
}
void DensePMF::prune (float thresh)
{
    int n=0;
    for (int i=0; i<MAX; ++i) {
        if (m_pmf[i] == 0) continue;
        if (m_pmf[i] < thresh) {
            m_pmf[i] = 0.0f;
            ++n;
        }
    }
    std::cout << n << " points pruned" << std::endl;
}
DensePMF product_and_sum (const DensePMF& other)
{//pointwise multiplication
    DensePMF result;
    DensePMF::const_iterator iter_lhs, iter_rhs, begin=other.begin();
    for (iter_lhs = begin; !iter_lhs.finished(); ++iter_lhs) {
        float p1 = iter_lhs.second();
        if (p1 == 0) continue;
        Point x1 = iter_lhs.first();
        for (iter_rhs = begin; !iter_rhs.finished(); ++iter_rhs) {
            float p2 = iter_rhs.second();
            if (p2 == 0) continue;
            Point x2 = iter_rhs.first();
    
            float p = p1 * p2;

            //product part
            Point x = x1 * x2;
            result[x] += p;

            //sum part
            x = x1 + x2;
            result[x] += p;
            
        }
    }
    return result;
}

int main (int argc, char** argv)
{
    typedef DensePMF PMF;

    float thresh = 1e-6; //default falue
    if (argc > 1) { thresh = atof(argv[1]); }

    std::cout << "grid size = " << X_MAX << " x " << Y_MAX << std::endl;
    std::cout << "threshold = " << thresh << std::endl;

    //int r0 = atoi(argv[2]);
    //int c0 = atoi(argv[3]);
    //Point x0(r0,c0);

    Point x0(X_POS, Y_POS);

    PMF points(x0);
    unsigned size=1, old_size=1;
    do {
        points *= 0.5;

        points = product_and_sum(points); //adds 1/4 tot for both product and sum
        points[x0] += 0.5;

        points.prune(thresh);

        old_size = size;
        size = points.size();
        std::cout << size << " points total" << std::endl;
        std::cout << "mass = " << points.mass() << std::endl;
    } while (old_size < size);

    points.writeToFile();

    return 0;
}

