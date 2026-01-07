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

def MakeIMatrix(cols, rows):
    """
    Function to create an identity matrix.
    :param cols: Number of columns.
    :param rows: Number of rows.
    :return: Identity matrix of size cols x rows.
    """
    return [[1 if x == y else 0 for y in range(cols)] for x in range(rows)]

def matrix_multiply(A, B):
    """
    Function that takes two matrices and multiplies them.
    :param A: First matrix.
    :param B: Second matrix.
    :return: The resulting matrix from multiplying A by B.
    """
    if len(A[0]) != len(B):
        raise ValueError("Matrix dimensions are incompatible for multiplication.")

    result = [[0 for _ in range(len(B[0]))] for _ in range(len(A))]

    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                result[i][j] += A[i][k] * B[k][j]

    return np.array(result)

def LU(A):
    """
    Function to perform LU decomposition.
    :param A: Matrix to be decomposed.
    :return: L and U matrices.
    """
    U = A
    L = MakeIMatrix(3, 3)
    for i in range(3):
        for j in range(3):
            if (i == 1 and j == 0 or i == 2 and j == 0 or i == 2 and j == 1):
                m_i = MakeIMatrix(3, 3)
                if L[j][j] != 0:
                    m_i[i][j] = -L[i][j] / L[j][j]

                U = matrix_multiply(m_i, U)
                L = matrix_multiply(L, get_inverse_matrix(m_i))
    return L, U

def get_inverse_matrix(A):
    """
    Function that takes a matrix and computes its inverse using elementary matrices.
    :param A: Matrix nxn.
    :return: Inverse of the matrix.
    """
    matrix_new_a = MakeIMatrix(3, 3)
    A_to_i = A
    for i in range(3):
        for j in range(3):
            if i == 1 and j == 0 or i == 2 and j == 0 or i == 2 and j == 1:
                m_i = MakeIMatrix(3, 3)
                if A_to_i[j][j] != 0:
                    m_i[i][j] = -A_to_i[i][j] / A_to_i[j][j]
                A_to_i = matrix_multiply(m_i, A_to_i)
                matrix_new_a = matrix_multiply(m_i, matrix_new_a)

    for i in range(2, -1, -1):
        for j in range(2, -1, -1):
            if i == 1 and j == 2 or i == 0 and j == 2 or i == 0 and j == 1:
                m_i = MakeIMatrix(3, 3)
                if A_to_i[j][j] != 0:
                    m_i[i][j] = -A_to_i[i][j] / A_to_i[j][j]
                A_to_i = matrix_multiply(m_i, A_to_i)
                matrix_new_a = matrix_multiply(m_i, matrix_new_a)

    for i in range(3):
        for j in range(3):
            if i == j:
                m_i = MakeIMatrix(3, 3)
                if A_to_i[i][j] != 0:
                    m_i[i][j] = 1 / A_to_i[i][j]
                A_to_i = matrix_multiply(m_i, A_to_i)
                matrix_new_a = matrix_multiply(m_i, matrix_new_a)
    return matrix_new_a

# הגדרת מטריצות A ו-b לצורך בדיקות
A = np.array([
    [2, 2, 2],
    [4, 2, 2],
    [2, 6, 2]
])
b = np.array([
    [1],
    [2],
    [3]
])

# הדפסת המטריצה ההפוכה של A
print_matrix(get_inverse_matrix(A))

# בדיקת תוצאה של כפל המטריצה ההפוכה של A ב-A
print_matrix(matrix_multiply(get_inverse_matrix(A), A))

# ביצוע פירוק LU על המטריצה A
l, u = LU(A)

# חישוב והדפסת המטריצות ההפוכות של L ו-U
l_1 = get_inverse_matrix(l)
u_1 = get_inverse_matrix(u)

# הדפסת תוצאה של כפל המטריצות ההפוכות של L ו-U ב-b
print_matrix(matrix_multiply(matrix_multiply(l_1, u_1), b))
