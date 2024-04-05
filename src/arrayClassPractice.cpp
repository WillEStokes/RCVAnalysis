#include <iostream>
#include <array>
#include <string.h>
// #include "arrayPractice.h"

template<typename T, size_t S>
class Array
{
public:
    constexpr int Size() const { return S; }

    T& operator[](size_t index) {return m_Data[index]; }
    const T& operator[](size_t index) const {return m_Data[index]; }

    T* Data() { return m_Data; }
    const T* Data() const { return m_Data; }
private:
    T m_Data[S];
};

// int arrayPractice(double* xdata, unsigned int element, double* vectElement)
int main()
{
    int size = 5;
    Array<int, 5> data;
    memset(&data[0], 0, data.Size() * sizeof(int));

    // const auto& arrayReference = data;
    for (int i = 0; i < data.Size() +1; i++)
    {
        // data[i] = 2;
        std::cout << data[i] << std::endl;
    }

    // printf("%s", "Hello World");
    std::cin.get();

    return 1;
}