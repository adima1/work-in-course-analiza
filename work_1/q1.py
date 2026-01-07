import subprocess
import sys

try:
    import numpy as np
except ImportError:
    print("numpy is not installed. Installing numpy...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "numpy"])
    import numpy as np

def print_matrix(matrix):
    """
    Function to print a matrix.
    :param matrix: Matrix to be printed.
    """
    for row in matrix:
        for element in row:
            print(element, end=" ")  # Print each element in the row
        print()  # Move to the next row
    print()


def matrix_multiply(A, B):
    """
    Function that takes two matrices and multiplies them
    :param A: First matrix
    :param B: Second matrix
    :return: The resulting matrix from multiplying A by B
    """
    if len(A[0]) != len(B):
        raise ValueError("Matrix dimensions are incompatible for multiplication.")

    result = [[0 for _ in range(len(B[0]))] for _ in range(len(A))]

    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                result[i][j] += A[i][k] * B[k][j]

    return np.array(result)


def MakeIMatrix(cols, rows):
    # Initialize a identity matrix
    return [[1 if x == y else 0 for y in range(cols)] for x in range(rows)]



def MaxNorm(matrix):
    """
    Function for calculating the max-norm of a matrix
    :param matrix: Matrix nxn
    :return:max-norm of a matrix
    """
    max_norm = 0
    for i in range(len(matrix)):
        norm = 0
        for j in range(len(matrix)):
            # Sum of organs per line with absolute value
            norm += abs(matrix[i][j])
        # Maximum row amount
        if norm > max_norm:
            max_norm = norm

    return max_norm

A = np.array([
    [1, -1, -2],
    [2, -3, -5],
    [-1, 3, 5]
])


def get_inverse_matrix(A):
    """
    Function that takes a matrix and computes its inverse using elementary matrices
    :param matrix: Matrix nxn
    :return:inverse of a matrix
    """
    matrix_new_a = MakeIMatrix(3, 3)
    A_to_i=A
    for i in range(3):
        for j in range(3):
            if(i==1 and j==0 or i==2 and j==0 or i==2 and j==1):
                m_i = MakeIMatrix(3, 3)
                m_i[i][j] = - A_to_i[i][j] /  A_to_i[j][j]
                A_to_i = matrix_multiply(m_i, A_to_i)
                matrix_new_a = matrix_multiply(m_i, matrix_new_a)

    for i in range(2,-1,-1):
        for j in range(2,-1,-1):
            if i==1 and j==2 or  i==0 and j==2 or  i==0 and j==1:
                m_i = MakeIMatrix(3, 3)
                m_i[i][j] = - A_to_i[i][j] /  A_to_i[j][j]
                A_to_i = matrix_multiply(m_i,A_to_i)
                matrix_new_a = matrix_multiply(m_i, matrix_new_a)

    for i in range(3):
        for j in range(3):
            if i==j:
                m_i = MakeIMatrix(3, 3)
                m_i[i][j]=1/ A_to_i[i][j]
                A_to_i=matrix_multiply(m_i, A_to_i)
                matrix_new_a=matrix_multiply(m_i,matrix_new_a)
    return matrix_new_a

def get_cond(A,B):
    """
    Function to calculate the COND of a matrix
    :param matrix: Matrix
    :return: The COND of the matrix
    """
    return MaxNorm(matrix_new_a)*MaxNorm(A)


matrix_new_a=get_inverse_matrix(A)
print_matrix(matrix_multiply(matrix_new_a,A))

print(get_cond(A,matrix_new_a))
