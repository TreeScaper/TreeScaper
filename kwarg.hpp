#pragma once
#define PRECISION double

class KwArg
{
private:
    // Users are responsible for type casting.
    // This class only keep pointers to existing data and never release their memory.
    // A copy of (the pointer of) data must be kept in order for releasing memories.
    unsigned int N;
    void **args;

public:
    KwArg() : N(0), args(nullptr){};
    KwArg(unsigned int len) : N(len), args(len ? new void *[len] : nullptr)
    {
        for (auto i = 0; i < N; i++)
            args[i] = nullptr;
    };
    KwArg(unsigned int len, void **KWARGS) : N(len), args(len ? new void *[len] : nullptr)
    {
        for (unsigned int i = 0; i < N; i++)
            args[i] = KWARGS[i];
    };
    KwArg(const KwArg &src) : KwArg(src.N, src.args){};
    ~KwArg()
    {
        // for (auto i = 0; i < N; i++)
        //     args[i] = nullptr;
        delete[] args;
    }

    friend void swap(KwArg &lhs, KwArg &rhs)
    {
        using std::swap;
        swap(rhs.N, lhs.N);
        swap(rhs.args, lhs.args);
    }

    KwArg &operator=(const KwArg &rhs)
    {
        KwArg temp(rhs);
        swap(*this, temp);
        return *this;
    }

    template <class T>
    void set(T *arg, unsigned int ind) { args[ind] = reinterpret_cast<void *>(arg); };

    void *operator[](int ind) { return args[ind]; };

    void *operator()(int ind) { return args[ind]; };
};

// void nullfunction(KwArg &kws, unsigned int *ind){};