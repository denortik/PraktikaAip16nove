#include <cstddef>
#include <iostream>
#include <stdexcept>

struct Mtx {
    size_t m;
    size_t n;
    double* data;

    Mtx()
        : m(0),
        n(0),
        data(nullptr) {
    }

    Mtx(size_t r, size_t c, double v = 0.0)
        : m(r),
        n(c),
        data(r* c != 0 ? new double[r * c] : nullptr) {
        for (size_t i = 0; i < r * c; ++i) {
            data[i] = v;
        }
    }

    ~Mtx() {
        delete[] data;
    }

    Mtx(const Mtx& o)
        : m(o.m),
        n(o.n),
        data(o.m* o.n != 0 ? new double[o.m * o.n] : nullptr) {
        for (size_t i = 0; i < m * n; ++i) {
            data[i] = o.data[i];
        }
    }

    Mtx& operator=(const Mtx& o) {
        if (this != &o) {
            size_t total = o.m * o.n;
            double* nw = total != 0 ? new double[total] : nullptr;
            for (size_t i = 0; i < total; ++i) {
                nw[i] = o.data[i];
            }
            delete[] data;
            data = nw;
            m = o.m;
            n = o.n;
        }
        return *this;
    }

    Mtx(Mtx&& o) noexcept
        : m(o.m),
        n(o.n),
        data(o.data) {
        o.m = 0;
        o.n = 0;
        o.data = nullptr;
    }

    Mtx& operator=(Mtx&& o) noexcept {
        steal(o);
        return *this;
    }

    void steal(Mtx& o) noexcept {
        if (this != &o) {
            delete[] data;
            m = o.m;
            n = o.n;
            data = o.data;
            o.m = 0;
            o.n = 0;
            o.data = nullptr;
        }
    }

    double& at(size_t i, size_t j) {
        return data[i * n + j];
    }

    double at(size_t i, size_t j) const {
        return data[i * n + j];
    }
};

Mtx add(Mtx a, Mtx b) {
    if (a.m != b.m || a.n != b.n) {
        throw std::invalid_argument("add: size mismatch");
    }
    Mtx r(a.m, a.n);
    for (size_t i = 0; i < r.m * r.n; ++i) {
        r.data[i] = a.data[i] + b.data[i];
    }
    return r;
}

Mtx mul(Mtx a, Mtx b) {
    if (a.n != b.m) {
        throw std::invalid_argument("mul: size mismatch");
    }
    Mtx r(a.m, b.n, 0.0);
    for (size_t i = 0; i < a.m; ++i) {
        for (size_t j = 0; j < b.n; ++j) {
            for (size_t p = 0; p < a.n; ++p) {
                r.at(i, j) += a.at(i, p) * b.at(p, j);
            }
        }
    }
    return r;
}

void print(const Mtx& mat) {
    for (size_t i = 0; i < mat.m; ++i) {
        for (size_t j = 0; j < mat.n; ++j) {
            std::cout << mat.at(i, j) << ' ';
        }
        std::cout << '\n';
    }
}

struct Mtxs {
    size_t k;
    Mtx* mtxs;

    Mtxs()
        : k(0),
        mtxs(nullptr) {
    }

    Mtxs(size_t n)
        : k(n),
        mtxs(n != 0 ? new Mtx[n] : nullptr) {
    }

    ~Mtxs() {
        delete[] mtxs;
    }

    Mtxs(const Mtxs& o)
        : k(o.k),
        mtxs(o.k != 0 ? new Mtx[o.k] : nullptr) {
        for (size_t i = 0; i < k; ++i) {
            mtxs[i] = o.mtxs[i];
        }
    }

    Mtxs& operator=(const Mtxs& o) {
        if (this != &o) {
            Mtxs tmp(o);
            swap(tmp);
        }
        return *this;
    }

    Mtxs(Mtxs&& o) noexcept
        : k(o.k),
        mtxs(o.mtxs) {
        o.k = 0;
        o.mtxs = nullptr;
    }

    Mtxs& operator=(Mtxs&& o) noexcept {
        if (this != &o) {
            delete[] mtxs;
            k = o.k;
            mtxs = o.mtxs;
            o.k = 0;
            o.mtxs = nullptr;
        }
        return *this;
    }

    void swap(Mtxs& o) noexcept {
        size_t tk = k;
        k = o.k;
        o.k = tk;

        Mtx* tm = mtxs;
        mtxs = o.mtxs;
        o.mtxs = tm;
    }

    size_t size() const {
        return k;
    }

    bool empty() const {
        return k == 0;
    }

    bool not_empty() const {
        return k != 0;
    }

    Mtx& operator[](size_t i) {
        return mtxs[i];
    }

    const Mtx& operator[](size_t i) const {
        return mtxs[i];
    }

    void erase(size_t from, size_t to) {
        if (from > to || to > k) {
            throw std::out_of_range("erase: bad range");
        }

        size_t cnt = to - from;
        if (cnt == 0) {
            return;
        }

        size_t newK = k - cnt;
        Mtx* nw = newK != 0 ? new Mtx[newK] : nullptr;

        for (size_t i = 0; i < from; ++i) {
            nw[i].steal(mtxs[i]);
        }
        for (size_t i = to; i < k; ++i) {
            nw[i - cnt].steal(mtxs[i]);
        }

        delete[] mtxs;
        mtxs = nw;
        k = newK;
    }

    void erase(size_t pos) {
        if (pos >= k) {
            throw std::out_of_range("erase: bad position");
        }
        erase(pos, pos + 1);
    }

    void erase_if(bool (*pred)(const Mtx&)) {
        size_t w = 0;
        for (size_t r = 0; r < k; ++r) {
            if (!pred(mtxs[r])) {
                if (w != r) {
                    mtxs[w].steal(mtxs[r]);
                }
                ++w;
            }
        }
        erase(w, k);
    }

    void insert(size_t pos, const Mtxs& other) {
        if (pos > k) {
            throw std::out_of_range("insert: bad position");
        }
        if (other.k == 0) {
            return;
        }
        if (&other == this) {
            Mtxs copy(other);
            insert(pos, copy);
            return;
        }

        Mtx* nw = new Mtx[k + other.k];

        for (size_t i = 0; i < pos; ++i) {
            nw[i].steal(mtxs[i]);
        }
        for (size_t i = 0; i < other.k; ++i) {
            nw[pos + i] = other.mtxs[i];
        }
        for (size_t i = pos; i < k; ++i) {
            nw[i + other.k].steal(mtxs[i]);
        }

        delete[] mtxs;
        mtxs = nw;
        k += other.k;
    }

    void insert(size_t pos, const Mtx& m) {
        Mtxs tmp(1);
        tmp.mtxs[0] = m;
        insert(pos, tmp);
    }

    void insert(size_t pos, size_t count, const Mtx& m) {
        Mtxs tmp(count);
        for (size_t i = 0; i < count; ++i) {
            tmp.mtxs[i] = m;
        }
        insert(pos, tmp);
    }

    Mtxs operator+(const Mtxs& rhs) const {
        Mtxs result(k + rhs.k);
        for (size_t i = 0; i < k; ++i) {
            result.mtxs[i] = mtxs[i];
        }
        for (size_t i = 0; i < rhs.k; ++i) {
            result.mtxs[k + i] = rhs.mtxs[i];
        }
        return result;
    }

    Mtx conv(Mtx m, Mtx(*t)(Mtx, Mtx)) const {
        for (size_t i = 0; i < k; ++i) {
            m = t(m, mtxs[i]);
        }
        return m;
    }
};

int main() {
    Mtx a(2, 2, 1.0);
    Mtx b(2, 2, 2.0);

    Mtxs arr;
    arr.insert(0, a);
    arr.insert(1, b);
    arr.erase(0);

    Mtxs sum = arr + arr;
    std::cout << "size: " << sum.size() << "\n";

    Mtx zero(2, 2, 0.0);
    print(sum.conv(zero, add));
    return 0;
}