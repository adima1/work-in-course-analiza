import numpy as np


def jacobi_method(coefficients, constants, initial_guess, tolerance=0.001, max_iterations=50):
    n = len(coefficients)
    x = initial_guess.copy()
    x_new = np.zeros(n)
    iteration = 0

    while iteration < max_iterations:
        for i in range(n):
            sum_ax = calculate_sum(coefficients, x, i)
            x_new[i] = (constants[i] - sum_ax) / coefficients[i, i]

        print(f"Jacobi Iteration {iteration + 1}: x = {x_new[0]:.6f}, y = {x_new[1]:.6f}, z = {x_new[2]:.6f}")

        if has_converged(x, x_new, tolerance):
            print(f'Jacobi method converged in {iteration + 1} iterations.')
            return x_new

        x = x_new.copy()
        iteration += 1

    print('Jacobi method did not converge within the specified tolerance and maximum iterations.')
    return x_new


def gauss_seidel_method(coefficients, constants, initial_guess, tolerance=0.001, max_iterations=50):
    n = len(coefficients)
    x = initial_guess.copy()
    iteration = 0

    while iteration < max_iterations:
        for i in range(n):
            sum_ax = calculate_sum(coefficients, x, i)
            x[i] = (constants[i] - sum_ax) / coefficients[i, i]

        print(f"Gauss-Seidel Iteration {iteration + 1}: x = {x[0]:.6f}, y = {x[1]:.6f}, z = {x[2]:.6f}")

        if has_converged_gauss_seidel(coefficients, constants, x, tolerance):
            print(f'Gauss-Seidel method converged in {iteration + 1} iterations.')
            return x

        iteration += 1

    print('Gauss-Seidel method did not converge within the specified tolerance and maximum iterations.')
    return x


def calculate_sum(coefficients, x, i):
    return np.dot(coefficients[i, :i], x[:i]) + np.dot(coefficients[i, i + 1:], x[i + 1:])


def has_converged(x, x_new, tolerance):
    return np.allclose(x, x_new, atol=tolerance)


def has_converged_gauss_seidel(coefficients, constants, x, tolerance):
    return np.linalg.norm(np.dot(coefficients, x) - constants) < tolerance


def is_diagonally_dominant(matrix):
    n = len(matrix)
    for i in range(n):
        row_sum = np.sum(np.abs(matrix[i])) - np.abs(matrix[i, i])
        if np.abs(matrix[i, i]) <= row_sum:
            return False
    return True


def solve_system(coefficients, constants, initial_guess):
    if is_diagonally_dominant(coefficients):
        print("Jacobi method:")
        jacobi_solution = jacobi_method(coefficients, constants, initial_guess)
        print("Jacobi Solution:", jacobi_solution)

        print("\nGauss-Seidel method:")
        gauss_seidel_solution = gauss_seidel_method(coefficients, constants, initial_guess)
        print("Gauss-Seidel Solution:", gauss_seidel_solution)
    else:
        print("Matrix is not diagonally dominant. Cannot guarantee convergence of Jacobi or Gauss-Seidel methods.")


def main():
    coefficients = np.array([[11, -5, 2],
                             [3, 12, -4],
                             [2, 1, 8]], dtype=float)
    constants = np.array([6, 25, -11], dtype=float)
    initial_guess = np.array([0, 0, 0], dtype=float)

    solve_system(coefficients, constants, initial_guess)


if __name__ == "__main__":
    main()
