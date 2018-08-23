This package implemented three algorithms for matrix transpose and multiplication, respectively: 1. brute force, 2. tiling only, 3. tiling plus multithreading.

To compile and generate the executable, type this in terminal window: g++ main.cpp -std=c++11 -lpthread
<Note that g++ main.cpp -std=c++11 -lthread does not compile in my mac nor ubuntu. lthread is not found>

This package is developed using Clion.  The project can be directly opened and run to obtain results comparing the algorithms.

The packaged is in my github repository: https://github.com/aureus5/MatrixMultiplicationTranspose.git  From here you can see my dev record.

commonUtils.h contains methods that are used by both matrix transpose and multiply.

Code for matrix transpose is in MatrixTranspose.h.

Code for Matrix multiply is in MatrixMultiplication.h.

For an example of how to use each algorithm for transpose, please refer to the driving method, void transposeMatrix() in MatrixTranspose.h. Parameters you can use include: the two dimensions of the source matrix, tile size, range of random numbers to fill the source matrix, and number of threads used in tile+thread method.

Likewise, for a reference to use the algorithms for multiplication, refer to the driver method void matrixMulplication() to see how to pass parameters such as the two dimensions of the two matrices, tile size, range of random numbers to fill the source matrix, and number of threads used in tile+thread method.

In the plotForComparingPerformance folder, I have plot the performance comparision for transpose and multiply using pyplot.  The source file is attached.  The plots are generated using python 3, in PyCharm.