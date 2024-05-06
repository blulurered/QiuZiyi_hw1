#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    int rows, cols; // the rows and columns of the matrix
    Matrix sum;     // the sum matrix of a and b
    int i, j;       // the iteration variables

    // check whether the size if incompatible
    if ((a.rows != b.rows) ||( a.cols != b.cols))
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    rows = a.rows;
    cols = a.cols;
    sum = create_matrix(rows, cols); // create a matrix
    for (i = 0; i < rows; i++)       // iterate over the rows
    {
        for (j = 0; j < cols; j++) // iterate over the columns
        {
            sum.data[i][j] = a.data[i][j] + b.data[i][j]; // calculate the sum for each element
        }
    }

    return sum;
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    int rows, cols; // the rows and columns of the matrix
    Matrix dif;     // the difference matrix of a and b
    int i, j;       // the iteration variables

    // check whether the size if incompatible
    if (a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    rows = a.rows;
    cols = a.cols;
    dif = create_matrix(rows, cols); // create a matrix
    for (i = 0; i < rows; i++)       // iterate over the rows
    {
        for (j = 0; j < cols; j++) // iterate over the columns
        {
            dif.data[i][j] = a.data[i][j] - b.data[i][j]; // calculate the sum for each element
        }
    }

    return dif;
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    int rows, cols; // the rows and columns of the matrix
    Matrix product; // the product matrix of a and b
    int i, j, k;    // the iteration variables

    // check whether the size if incompatible
    if (a.cols != b.rows)
    {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    rows = a.rows;
    cols = b.cols;
    product = create_matrix(rows, cols); // create a matrix
    for (i = 0; i < rows; i++)           // iterate over the rows
    {
        for (j = 0; j < cols; j++) // iterate over the columns
        {
            product.data[i][j] = 0; // initialize the elements to 0
            for (k = 0; k < a.cols; k++)
                product.data[i][j] += a.data[i][k] * b.data[k][j]; // calculate the product for each element
        }
    }

    return product;
}

Matrix scale_matrix(Matrix a, double k)
{
    int rows, cols; // the rows and columns of the matrix
    Matrix product; // the product matrix of a and double k
    int i, j;       // the iteration variables

    rows = a.rows;
    cols = a.cols;
    product = create_matrix(rows, cols); // create a matrix
    for (i = 0; i < rows; i++)           // iterate over the rows
    {
        for (j = 0; j < cols; j++) // iterate over the columns
        {
            product.data[i][j] = a.data[i][j] * k; // calculate the product for each element
        }
    }

    return product;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix transposed;

    // exchange the rows and columns
    transposed.rows = a.cols;
    transposed.cols = a.rows;
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            transposed.data[j][i] = a.data[i][j]; // calculate the element
        }
    }

    return transposed;
}

// the main idea is to design a recursive algorithm
// we calculate the determinant by expanding the first column
double det_matrix(Matrix a)
{
    double determinant = 0.0;
    int sign = 1;
    Matrix submatrix;
    int i, j, k;

    // if a is not a square matrix
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    // the basic case, a is a 1*1 matrix
    if (a.rows == 1)
    {
        return a.data[0][0];
    }
    submatrix.rows = a.rows - 1;
    submatrix.cols = a.cols - 1;
    for (i = 0; i < a.rows; i++)
    {
        // for each element in the first column, caculate the corresponding submatrix
        for (j = 0; j < a.rows - 1; j++)
        {
            if (j < i)
            {
                for (k = 0; k < a.cols - 1; k++)
                    submatrix.data[j][k] = a.data[j][k + 1];
            }
            else
            {
                for (k = 0; k < a.cols - 1; k++)
                    submatrix.data[j][k] = a.data[j + 1][k + 1];
            }
        }
        // add the product of element and the determinant of its submatrix
        // the determinant of its submatrix is calculated recursively
        determinant += sign * a.data[i][0] * det_matrix(submatrix);
        sign = -sign;
    }

    return determinant;
}

// the main idea is calculating through the Adjugate of matrix a
Matrix inv_matrix(Matrix a)
{
    Matrix inverse;
    Matrix submatrix;                   // the submatrix of each element in a
    double determinant = det_matrix(a); // the determinant of a
    int i, j;
    int row, col;

    // check whether a is a square matrix
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    // check whether a is singular
    if (determinant == 0)
    {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    inverse = create_matrix(a.rows, a.cols);
    submatrix = create_matrix(a.rows - 1, a.cols - 1);
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            int subrow = 0, subcol = 0;
            for (row = 0; row < a.rows; row++)
            {
                if (row == i) // skip a row if
                    continue;
                for (col = 0; col < a.cols; col++)
                {
                    if (col == j)
                        continue;
                    submatrix.data[subrow][subcol] = a.data[row][col]; // calculate the element according to a
                    subcol++;
                }
                subrow++;
                subcol = 0; // reset the subcol
            }
            double cofactor = pow(-1, i + j) * det_matrix(submatrix); // calculate the cofactors
            inverse.data[i][j] = cofactor / determinant;
        }
    }
    inverse = transpose_matrix(inverse); // transpose the matrix

    return inverse;
}

// the main idea is to perform the Gaussian elimination
int rank_matrix(Matrix a)
{
    int rank = 0;
    int i, j, k, col;
    Matrix temp = a; // we will operate on the matrix temp, a remains unchanged

    for (j = 0; j < a.cols; j++)
    {                                                            // iterate over each columns
        int pivotRow = rank;                                     // initialize the pivotRow
        while (pivotRow < a.rows && temp.data[pivotRow][j] == 0) // iterate until finding a nonzero pivot
            pivotRow++;
        if (pivotRow == a.rows) // the case when all the elements are zero
            continue;
        // Swap the row of rank and pivotRow
        for (int col = 0; col < a.cols; col++)
        {
            double tempVal = temp.data[rank][col];
            temp.data[rank][col] = temp.data[pivotRow][col];
            temp.data[pivotRow][col] = tempVal;
        }
        // Perform row operations to make all elements below the pivot element zero
        for (int i = rank + 1; i < a.rows; i++)
        {
            double factor = temp.data[i][j] / temp.data[rank][j];
            for (k = j; k < a.cols; k++)
            {
                temp.data[i][k] -= factor * temp.data[rank][k];
            }
        }
        rank++;
    }

    return rank;
}

double trace_matrix(Matrix a)
{
    // if a is not a square matrix
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }

    double trace = 0.0;
    // add the elements on the diagonal
    for (int i = 0; i < a.rows; i++)
    {
        trace += a.data[i][i];
    }

    return trace;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}