#include <QCoreApplication>

#include "Math/matrix.h"
#include <iostream>
int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    double arr[2][2] = {{1,2}, {3,4}};
    double arr2[2][2] = {{1,1},{1,1}};
    matrix m_arr(arr, 2, 2);
    matrix m_arr2(arr2, 2, 2);

    std::cout << m_arr;
    return a.exec();
}
