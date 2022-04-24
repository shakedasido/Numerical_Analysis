def Jacobi(mtx, vector_b):
    """
        Jacobi Method for solving a converge matrix.

        Parameters
        ----------
        mtx: nested list (nXn)
        vector_b: nested list (nX1)

        Returns
        -------
        None

    """
    x_start, y_start, z_start, count, epsilon, iterate_if_not_dom = 0, 0, 0, 0, 0.000001, 100
    flag = True
    check = CheckDominantDiag(mtx)
    if not check:
        mtx, vector_b = pivoting(mtx, vector_b)
        if not CheckDominantDiag(mtx):
            print("The matrix does not have dominant diagonal. Let's see if it converges")
            flag = False
    else:
        print("The matrix have dominant diagonal")

    while True:
        count += 1
        x = (vector_b[0][0] - mtx[0][1] * y_start - mtx[0][2] * z_start) / mtx[0][0]
        y = (vector_b[1][0] - mtx[1][0] * x_start - mtx[1][2] * z_start) / mtx[1][1]
        z = (vector_b[2][0] - mtx[2][0] * x_start - mtx[2][1] * y_start) / mtx[2][2]
        print(f'#{count} iteration: x = {x}\ty = {y}\tz = {z}')

        # if the matrix does not have dominant diagonal, and it didn't converge after 100 iterations, then:
        if count == iterate_if_not_dom and not flag:
            print("The system does not converge!\n")
            return

        if abs(x - x_start) > epsilon:
            x_start = x
            y_start = y
            z_start = z
        else:
            break
    print("")
    if not flag:
        print("Although the matrix does not have dominant diagonal, ")
    print(f'Result are: x = {x}\ty = {y}\tz = {z}\n')


def GaussSeidel(mtx, vector_b):
    """
        Gauss Seidel Method for solving a converge matrix.

        Parameters
        ----------
        mtx: nested list (nXn)
        vector_b: nested list (nX1)

        Returns
        -------
        None

    """
    x_start, y_start, z_start, count, epsilon, iterate_if_not_dom = 0, 0, 0, 0, 0.000001, 100
    flag = True
    check = CheckDominantDiag(mtx)
    if not check:
        mtx, vector_b = pivoting(mtx, vector_b)
        if not CheckDominantDiag(mtx):
            print("The matrix does not have dominant diagonal. Let's see if it converges")
            flag = False
    else:
        print("The matrix have dominant diagonal")

    while True:
        count += 1
        x = (vector_b[0][0] - mtx[0][1] * y_start - mtx[0][2] * z_start) / mtx[0][0]
        y = (vector_b[1][0] - mtx[1][0] * x - mtx[1][2] * z_start) / mtx[1][1]
        z = (vector_b[2][0] - mtx[2][0] * x - mtx[2][1] * y) / mtx[2][2]
        print(f'#{count} iteration: x = {x}\ty = {y}\tz = {z}')

        # if the matrix does not have dominant diagonal, and it didn't converge after 100 iterations, then:
        if count > iterate_if_not_dom and not flag:
            print("The system does not converge!\n")
            return

        if abs(x - x_start) > epsilon:
            x_start = x
            y_start = y
            z_start = z
        else:
            break
    print("")
    if not flag:
        print("Although the matrix does not have dominant diagonal, ")
    print(f'Result are: x = {x}\ty = {y}\tz = {z}\n')


def CheckDominantDiag(matrix):
    """
        A function that checks if a matrix has a dominant diagonal.

        Parameters
        ----------
        matrix: nested list (nXn)

        Returns
        -------
        Boolean Result

    """
    for i in range(len(matrix)):
        Sum = 0
        for j in range(len(matrix)):
            if matrix[i][j] != (matrix[i][i]):
                Sum += abs(matrix[i][j])
        if abs(matrix[i][i]) < Sum:
            return False
    return True


def print_matrix(mtx):
    """
        gets a matrix and print it

        Parameters
        ----------
        mtx : nested list (nXn)

        Returns
        -------
        None
    """
    for line in mtx:
        print("[" + '  '.join(map(str, line)) + "]")


def mul_mtx_with_mtx(mtx1, mtx2):
    """
        multiple one matrix in another

        Parameters
        ----------
        mtx1: nested list (nXn) (elementary)
        mtx2: nested list (nXn)

        Returns
        -------
        The new matrix out of the multiplication.

    """
    new_matrix = []
    a = 0
    for i in range(len(mtx1)):
        new_matrix_row = []
        for j in range(len(mtx2[0])):
            for k in range(len(mtx2)):
                a += mtx1[i][k] * mtx2[k][j]
            new_matrix_row.append(a)
            a = 0
        new_matrix.append(new_matrix_row)
    return new_matrix


def create_mtx_i(rows, cols):
    """
        create a new identity matrix

        Parameters
        ----------
        rows: int
        cols: int

        Returns
        -------
        The new identity matrix.

    """

    i_mtx = []
    for i in range(rows):
        new_matrix_row = []
        for j in range(cols):
            if i == j:
                new_matrix_row.append(1)
            else:
                new_matrix_row.append(0)
        i_mtx.append(new_matrix_row)
    return i_mtx


def pivoting(mtx, vector_b):
    """
        func that replace rows.

        Parameters
        ----------
        mtx : nested list (nXn)
        vector_b: list

        Returns
        -------
        matrix after changing rows pivoting process.

    """
    for i in range(len(mtx)):
        pivot = abs(mtx[i][i])
        for j in range(i, len(mtx)):
            if abs(mtx[j][i]) > pivot:
                elementary_switched = create_mtx_i(len(mtx[0]), len(mtx))
                # switch the rows in elementary_switched matrix in order to have bigger pivot
                temp_row_el = elementary_switched[j]
                elementary_switched[j] = elementary_switched[i]
                elementary_switched[i] = temp_row_el
                # mul the matrix with the elementary matrix and print it
                print("\nAfter the pivoting process\n")
                mtx = mul_mtx_with_mtx(elementary_switched, mtx)
                vector_b = mul_mtx_with_mtx(elementary_switched, vector_b)
                pivot = abs(mtx[i][i])

    print("After changing the rows:\nMatrix: \n")
    print_matrix(mtx)
    print("Vector B:\n")
    print_matrix(vector_b)
    return [mtx, vector_b]


print("\nMatrix:\n")
matrixA = [[4, 2, 0], [2, 10, 4], [0, 4, 5]]
print_matrix(matrixA)
print("\nVector B:\n")
vectorB = [[2], [6], [5]]
print_matrix(vectorB)
print("\n")
while True:
    x = int(input('Please choose the method:\n1.Jacobi\n2.GaussSeidel\n'))
    if x == 1 or x == 2:
        break
    print("try again")

if x == 1:
    print("Solving by Jacobi method:")
    Jacobi(matrixA, vectorB)
else:
    print("Solving by GaussSeidel method:")
    GaussSeidel(matrixA, vectorB)





