import sys
import subprocess

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Check if numpy is installed
try:
    import numpy as np
except ImportError:
    print("Numpy not found, installing...")
    install('numpy')
    import numpy as np

# Define bcolors class for terminal colors
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#-------------------------------------------------------->polynomial
def swap_row(mat, i, j):
    for k in range(len(mat[0])):
        mat[i][k], mat[j][k] = mat[j][k], mat[i][k]

def forward_substitution(mat):
    N = len(mat)
    for k in range(N):

        # Partial Pivoting: Find the pivot row with the largest absolute value in the current column
        pivot_row = k
        v_max = abs(mat[pivot_row][k])
        for i in range(k + 1, N):
            if abs(mat[i][k]) > v_max:
                v_max = abs(mat[i][k])
                pivot_row = i

        # If the principal diagonal element is zero, the matrix is singular
        if mat[pivot_row][k] == 0:
            return k  # Matrix is singular

        # Swap the current row with the pivot row
        if pivot_row != k:
            swap_row(mat, k, pivot_row)

        for i in range(k + 1, N):
            # Compute the multiplier
            m = mat[i][k] / mat[k][k]

            # Subtract the multiple of the corresponding row element
            for j in range(k + 1, N + 1):
                mat[i][j] -= mat[k][j] * m

            # Set the lower triangular matrix element to zero
            mat[i][k] = 0

    return -1

def backward_substitution(mat):
    N = len(mat)
    x = np.zeros(N)  # An array to store the solution

    # Start calculating from the last equation up to the first
    for i in range(N - 1, -1, -1):
        x[i] = mat[i][N]

        # Initialize j to i+1 since matrix is upper triangular
        for j in range(i + 1, N):
            x[i] -= mat[i][j] * x[j]

        x[i] = x[i] / mat[i][i]

    return x

def gaussianElimination(mat):
    singular_flag = forward_substitution(mat)

    if singular_flag != -1:
        if mat[singular_flag][-1] != 0:
            return "Singular Matrix (Inconsistent System)"
        else:
            return "Singular Matrix (May have infinitely many solutions)"

    # If matrix is non-singular: get solution to system using backward substitution
    return backward_substitution(mat)

def polynomialInterpolation(matrix, table_points, x):
    b = [[point[1]] for point in table_points]
    matrix_new = np.hstack((matrix, b))
    matrix_sol = gaussianElimination(matrix_new)

    if isinstance(matrix_sol, str):
        print(matrix_sol)
        return None

    result = sum([matrix_sol[i] * (x ** i) for i in range(len(matrix_sol))])
    print(bcolors.OKBLUE, "\nThe polynomial:", bcolors.ENDC)
    print('P(X) = ' + ' + '.join([f'({matrix_sol[i]}) * x^{i}' for i in range(len(matrix_sol))]))
    print(bcolors.OKGREEN, f"\nThe Result of P(X={x}) is:", bcolors.ENDC, result)

    return result

def Prerequisite(table_points):
    matrix = [[point[0] ** i for i in range(len(table_points))] for point in table_points]  # Makes the initial matrix

    if np.isclose(np.linalg.det(matrix), 0):
        print("Singular Matrix")
        return None

    return matrix

#------------------------------------------------------------>linear

def linearInterpolation(table_points, point):
    length = len(table_points)
    for i in range(length - 1):
        if table_points[i][0] <= point <= table_points[i + 1][0]:
            x1, y1 = table_points[i]
            x2, y2 = table_points[i + 1]
            result = (((y1 - y2) / (x1 - x2)) * point) + ((y2 * x1) - (y1 * x2)) / (x1 - x2)
            print(bcolors.BOLD, "\nThe approximation (interpolation) of the point", point, "is:", bcolors.ENDC, round(result, 4))
            return
    x1, y1 = table_points[-2]
    x2, y2 = table_points[-1]
    m = (y1 - y2) / (x1 - x2)
    b = y1 - m * x1
    result = m * point + b
    print(bcolors.BOLD, "\nThe approximation (extrapolation) of the point", point, "is:", bcolors.ENDC, round(result, 4))

#----------------------------------------------------------->lagrange

def lagrange_interpolation(x_data, y_data, x):
    """
    Lagrange Interpolation

    Parameters:
    x_data (list): List of x-values for data points.
    y_data (list): List of y-values for data points.
    x (float): The x-value where you want to evaluate the interpolated polynomial.

    Returns:
    float: The interpolated y-value at the given x.
    """
    n = len(x_data)
    result = 0.0
    expression = ""

    for i in range(n):
        term = y_data[i]
        expression_term = f"({y_data[i]})"
        for j in range(n):
            if i != j:
                term *= (x - x_data[j]) / (x_data[i] - x_data[j])
                expression_term += f" * ((x - {x_data[j]}) / ({x_data[i]} - {x_data[j]}))"
        result += term
        expression += f" + {expression_term}" if i > 0 else expression_term

    print(bcolors.OKGREEN, "\nInterpolating polynomial P(x) =", expression, bcolors.ENDC)

    return result, expression

#----------------------------------------------------------->main

if __name__ == '__main__':
    table_points = [(1, 3), (2, 4), (3, -1)]
    x_data = [1, 2, 3]
    y_data = [3, 4, -1]
    x = 1.5
    choice = int(input("Please enter: \n1- Linear Interpolation \n2- Polynomial Interpolation \n3- Lagrange Interpolation\n"))

    if choice == 1:
        print(bcolors.OKBLUE, "----------------- Interpolation & Extrapolation Methods -----------------\n", bcolors.ENDC)
        print(bcolors.OKBLUE, "Table Points: ", bcolors.ENDC, table_points)
        print(bcolors.OKBLUE, "Finding an approximation to the point: ", bcolors.ENDC, x)
        linearInterpolation(table_points, x)
        print(bcolors.OKBLUE, "\n---------------------------------------------------------------------------\n", bcolors.ENDC)

    elif choice == 2:
        print(bcolors.OKBLUE, "----------------- Interpolation & Extrapolation Methods -----------------", bcolors.ENDC)
        matrix = Prerequisite(table_points)
        if matrix is not None:
            print(bcolors.OKBLUE, "Table Points: ", bcolors.ENDC, table_points)
            print(bcolors.OKBLUE, "Finding an approximation to the point: ", bcolors.ENDC, x)
            polynomialInterpolation(matrix, table_points, x)
        print(bcolors.OKBLUE, "---------------------------------------------------------------------------", bcolors.ENDC)

    elif choice == 3:
        y_interpolate, p = lagrange_interpolation(x_data, y_data, x)
        print(bcolors.OKBLUE, "\nInterpolated value at x =", x, "is y =", y_interpolate, bcolors.ENDC)
