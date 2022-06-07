import math
from tabulate import tabulate


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


def get_matrix_inverse(mtx, vector_b):
    """
        function that gets matrix and inverse it

        Parameters
        ----------
        mtx: nested list (nXn) (elementary)
        vector_b: list

        Returns
        -------
        The inverse matrix.

    """
    # we assume that the matrix is invertible
    # result matrix initialized as singularity matrix

    result = create_mtx_i(len(mtx), len(mtx))

    # row by row
    for i in range(len(mtx[0])):
        # pivoting: check if the pivot is 0/bigger number then one/1 and make changes, in order to make pivot=1
        mtx, vector_b = switch_rows(mtx, vector_b)
        elementary = create_mtx_i(len(mtx[0]), len(mtx))
        elementary[i][i] = 1 / mtx[i][i]
        # gets the result out of the multiplication of the elementary matrix with the last result.
        result = mul_mtx_with_mtx(elementary, result)
        mtx = mul_mtx_with_mtx(elementary, mtx)
        # turn to zero all the elements below the pivot
        for j in range(i + 1, len(mtx)):
            elementary = create_mtx_i(len(mtx[0]), len(mtx))
            elementary[j][i] = -(mtx[j][i])
            mtx = mul_mtx_with_mtx(elementary, mtx)
            result = mul_mtx_with_mtx(elementary, result)
    # turn to zero all the elements above the pivot
    for i in range(len(mtx[0]) - 1, 0, -1):
        for j in range(i - 1, -1, -1):
            elementary = create_mtx_i(len(mtx[0]), len(mtx))
            elementary[j][i] = -(mtx[j][i])
            mtx = mul_mtx_with_mtx(elementary, mtx)
            result = mul_mtx_with_mtx(elementary, result)

    return result


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


def switch_rows(mtx, vector_b):
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
                # mul the matrix with the elementary matrix
                mtx = mul_mtx_with_mtx(elementary_switched, mtx)
                vector_b = mul_mtx_with_mtx(elementary_switched, vector_b)
                pivot = abs(mtx[i][i])

    return [mtx, vector_b]


def switch_rows_if_zero(mtx, vector_b):
    """
        func that replace rows if one of the rows is zero.

        Parameters
        ----------
        mtx: nested list (nXn)
        vector_b: list

        Returns
        -------
        matrix after changing rows pivoting process.

    """
    for i in range(len(mtx)):
        for j in range(i, len(mtx)):
            # switch the rows in matrix in order to have non zero pivot
            temp_row_mtx = mtx[j]
            mtx[j] = mtx[i]
            mtx[i] = temp_row_mtx
            # switch the rows in the b vector if we have zero pivot.
            temp_b = vector_b[j]
            vector_b[j] = vector_b[i]
            vector_b[i] = temp_b

    return [mtx, vector_b]


def u_and_l_matrix(mtx, vector_b):
    """
        function that convert the matrix into two matrix's:  Upper triple matrix and Lower triple matrix

        Parameters
        ----------
        mtx: nested list (nXn)
        vector_b: list

        Returns
        -------
        Upper triple matrix

    """

    # let's switch rows in order to get the bigger pivot
    mtx, vector_b = switch_rows(mtx, vector_b)

    l_mtx = create_mtx_i(len(mtx[0]), len(mtx))
    # row by row
    for i in range(len(mtx[0])):
        # pivoting: check if the pivot is 0/bigger number then one/1 and make changes, in order to make pivot=1
        if mtx[i][i] == 0:
            mtx, vector_b = switch_rows_if_zero(mtx, vector_b)
        for j in range(i + 1, len(mtx)):
            if mtx[j][i] != 0 and mtx[j][i] != 0.0 and mtx[j][i] != -0.0:
                # create another i matrix
                elementary = create_mtx_i(len(mtx[0]), len(mtx))
                # make's all the elements under the pivot- 0
                # change the element in the i matrix in order to create a proper elementary matrix
                elementary[j][i] = -(mtx[j][i]) / mtx[i][i]
                l_mtx[j][i] = -(elementary[j][i])
                mtx = mul_mtx_with_mtx(elementary, mtx)

    return [mtx, l_mtx]


def solve_matrix(mtx, vector_b):
    """
         function that solve the matrix

        Parameters
        ----------
        mtx: nested list (nXn)
        vector_b: list

        Returns
        -------
        the result of the matrix - x

        """
    u_mtx, l_mtx = u_and_l_matrix(mtx, vector_b)

    # solving the matrix and get x, the result.
    # U(-)*L(-1)*VECTOR_B = X:

    U_inverse = get_matrix_inverse(u_mtx, vector_b)
    L_inverse = get_matrix_inverse(l_mtx, vector_b)

    # Create U^(-1)*l^(-1)*VECTOR_B = x:
    x = mul_mtx_with_mtx(U_inverse, mul_mtx_with_mtx(L_inverse, vector_b))

    return x


def create_mtx_spline(size, lamda, mu):
    """
        create a new matrix

        Parameters
        ----------
        size: int
        lamda: list
        mu: list

        Returns
        -------
        The new matrix for the spline interpolation.


    """

    i_mtx = []

    if not lamda and not mu:
        for i in range(size):
            new_matrix_row = []
            for j in range(size):
                if i == j:
                    new_matrix_row.append(2)
                elif abs(i - j) == 1:
                    new_matrix_row.append(1)
                else:
                    new_matrix_row.append(0)
            i_mtx.append(new_matrix_row)
    else:
        for i in range(size):
            new_matrix_row = []
            for j in range(size):
                if i == j:
                    new_matrix_row.append(2)
                elif (i - j) == 1:
                    new_matrix_row.append(mu[j])
                elif (j - i) == 1:
                    new_matrix_row.append(lamda[i])
                else:
                    new_matrix_row.append(0)
            i_mtx.append(new_matrix_row)
    return i_mtx


def linear_interpolation(table, x):
    """
    A linear interpolation connects a straight line between any two consecutive information points
    and in order to know the value of the new point we will use the straight line equation

        Parameters
        ----------
        table : Points table (a straight line will be connected between them)
        x : The x of the point we want to find its f(X)

        Returns
        -------
        None

    """
    name, flag, f_of_x = "Linear", False, 0
    Print_table(table)

    for i in range(len(table) - 1):
        if table[i][0] <= x <= table[i + 1][0] or table[i][0] > x > table[i + 1][0]:
            f_of_x = (((table[i][1] - table[i + 1][1]) / (table[i][0] - table[i + 1][0])) * x) + (
                    ((table[i + 1][1] * table[i][0]) - (table[i][1] * table[i + 1][0])) / (
                        table[i][0] - table[i + 1][0]))
            print(f"The result of the {name} interpolation on {x} is: f({x}) = {round(f_of_x, 5)}")
            flag = True
    if not flag:
        print("There are no values that satisfy a <= x <= b or a > x > b")


def polynomial_interpolation(table, x):
    """
    A polynomial interpolation- If n + 1 points are given, we can construct a polynomial from above n that passes through these points

    Parameters
    ----------
    table : Points table
    x : The x of the point we want to find its f(X)

    Returns
    -------
    None

    """

    name = "Polynomial"
    Print_table(table)
    # initial matrix
    matrix = [[point[0] ** i for i in range(len(table))] for point in table]
    # vector b
    vector_b = [[point[1]] for point in table]

    A = solve_matrix(matrix, vector_b)

    p_of_x = sum([A[i][0] * (x ** i) for i in range(len(A))])
    print(f"The result of the {name} interpolation on {x} by polynomial of {len(A)-1} is: p({x}) = {p_of_x}")


def neville_interpolation(table, x):
    """
    Neville interpolation -We will take 4 points from the table that are close to the point we get.
    This algorithm relies on solutions of induction polynomials.

        Parameters
        ----------
        table : table of points
        x : The x of the point we want to find its f(X)

        Returns
        -------
        None
    """
    name = "Neville"
    Print_table(table)
    for i in range(len(table)):
        if x == table[i][0]:
            print("The value of ", x, " is: ", round(table[i][1], 6))
            return

    tempTable = [[], [], [], []]
    k = 0
    for i in range(len(table) - 1):
        if table[i][0] <= x <= table[i + 1][0]:
            print(
                "The new point is between the values: " + str(table[i][0]) + " to " + str(table[i + 1][0]))
            k = i - 1
            if k<0:
                k=0
    for j in range(4):
        tempTable[j].append(table[k][0])
        tempTable[j].append(table[k][1])
        k = k + 1

    P01 = solve_neville(x, tempTable[0][0], tempTable[1][0], tempTable[0][1], tempTable[1][1])
    P12 = solve_neville(x, tempTable[1][0], tempTable[2][0], tempTable[1][1], tempTable[2][1])
    P23 = solve_neville(x, tempTable[2][0], tempTable[3][0], tempTable[2][1], tempTable[3][1])

    P02 = solve_neville(x, tempTable[0][0], tempTable[2][0], P01, P12)
    P13 = solve_neville(x, tempTable[1][0], tempTable[3][0], P12, P23)
    P03 = solve_neville(x, tempTable[0][0], tempTable[3][0], P02, P13)

    print(f"The result of the {name} interpolation on {x}  is: p({x}) = {round(P03, 5)}")


def solve_neville(x, x1, x2, p1, p2):
    """
    A Neville interpolation- the equation

        Parameters
        ----------
        x : The x of the point we want to find its f(X)
        x1 : x of point from the table
        x2 : x of point from the table
        p1 : the result of polynomials 1
        p2 : the result of polynomials 2

        Returns
        -------
        result of the equation

    """
    result = 0
    result = ((x - x1) * p2 - (x - x2) * p1) / (x2 - x1)
    return result


def spline_interpolation(table, x):
    """
    spline interpolation - in given table and x value, and defined elements for natural spline/
    regular spline. we create a matrix, calculate values and equations and finds f(x).

    Parameters
    ----------
    table : table of points
    x : The x of the point we want to find its f(X)

    Returns
    -------
    None
    """
    h_of_i = []
    # The distance between 2 x's : h(i) = x(i+1)-x(i), for i = 0, .., n-1
    for i in range(len(table) - 1):
        h_of_i.append(table[i + 1][0] - table[i][0])

    Print_Spline_table(table)
    choice = input("\nPlease choose:\n"
                   "1.Natural Cubic Spline\n"
                   "2.Cubic Spline\n")
    if choice == "1":
        d_0, d_n, mu_n, lam_0, name = 0, 0, 0, 0, "Natural Cubic Spline"
    else:
        choice = input("Please choose:\n"
                       "1.Get default f'(0) and f'(n)\n"
                       "2.Insert self parameters\n")
        if choice == "2":
            f_0 = float(input("f'(0): "))
            f_n = float(input("f'(0): "))
        else:
            f_0 = 1
            f_n = 0
        mu_n, lam_0, name = 1, 1, "Cubic Spline"
        d_0 = ((6 / h_of_i[0]) * (((table[1][1] - table[0][1]) / h_of_i[0]) - f_0))
        d_n = ((6 / h_of_i[len(table) - 2]) * (
                f_n - ((table[len(table) - 1][1] - table[len(table) - 2][1]) / h_of_i[0])))

    solve_Spline(table, x, d_0, d_n, mu_n, lam_0, h_of_i, name)


def solve_Spline(table, x, d_0, d_n, mu_n, lam_0, h_of_i, name):
    """
        A spline interpolation- the creation of the matrix and the solve of f(x)

            Parameters
            ----------
            table : table of points
            x : The x of the point we want to find its f(X)
            d_0 : int / float
            d_n : int / float
            mu_n : int / float
            lam_0 : 0/1
            h_of_i : int / float
            name : string

            Returns
            -------
            result of the equation

        """
    mtx = create_mtx_spline(len(table), [], [])
    lamda, mu, d = [], [], []
    # i made mu(1) to be mu(0) because the range of a list starts from 0 and not 1
    # mu(n)=mu(n-1) now. so lamda(0)=mu(n-1)=0 in natural spline, and lamda(0)=mu(n-1)=1 in regular spline.
    lamda.append(lam_0), d.append([round(d_0, 6)])

    # lamda(a) = h_of_i(a) / (h_of_i(a) + h_of_i(a-1), mu(a) = 1 - lamda(a) , for all a=1,..., n-1
    # mu(1,...,n), lamda(0,...,n-1)--> for the project: mu(0,...,n-1), lamda(0,...,n-1) --> mu(a) = 1-lamda(a+1)
    for i in range(1, len(mtx) - 1):  # len(mtx) - 1 isn't include, so until len(mtx) - 2
        lamda.append(h_of_i[i] / (h_of_i[i] + h_of_i[i - 1]))
        d.append([round((6 / (h_of_i[i] + h_of_i[i - 1])) * (((table[i + 1][1] - table[i][1]) / h_of_i[i]) -
                                                             ((table[i][1] - table[i - 1][1]) / h_of_i[i - 1])), 6)])
        if i <= len(mtx) - 2:
            mu.append(1 - lamda[i])

    mu.append(mu_n), d.append([round(d_n, 6)])

    for i in range(len(table)):
        for j in range(len(table)):
            if (i - j) == 1:
                mtx[i][j] = round(mu[j], 10)
            elif (j - i) == 1:
                mtx[i][j] = round(lamda[i], 10)

    M = solve_matrix(mtx, d)
    for i in range(len(M)):
        if M[i] == [0] or M[i] == [0.0] or M[i] == [-0.0]:
            M[i] = [0]
        else:
            M[i] = [round(M[i][0], 6)]


    # if the list is empty or The system does not converge.
    if not M:
        print("Unable to solve the matrix. No Available interpolation foe that table.")

    s_of_x = 0
    for i in range(len(table) - 1):
        if table[i][0] < x < table[i + 1][0] or table[i][0] > x > table[i + 1][0]:
            a = table[i][0]
            b = table[i + 1][0]
            s_of_x = ((((pow(b - x, 3) * M[i][0]) + (pow(x - a, 3) * M[i + 1][0])) / (6 * h_of_i[i])) + (
                    (((b - x) * table[i][1]) + ((x - a) * table[i + 1][1])) / h_of_i[i]) - (
                              ((b - x) * M[i][0]) + ((x - a) * M[i + 1][0])) * (h_of_i[i] / 6))
            break

    print(f"The result of the {name} interpolation on {x} is: {s_of_x}")


def lagrange_interpolation(table, x):
    """
    lagrange interpolation -look for a polynomial from above n that passes through some point and resets at all other points
    A polynomial is found by connecting polynomials and not by a matrix.
    This algorithm relies on solutions of induction polynomials.

    Parameters
    ----------
    table : table of points
    x : The x of the point we want to find its f(X)

    Returns
    -------
    None
    """

    name = "Lagrange"
    Print_table(table)
    x_list, y_list = [], []
    p_of_n = 0
    # collect all the x's from the table as a list
    # collect all the y's from the table as a list
    for i in range(len(table)):
        x_list.append(table[i][0])
        y_list.append(table[i][1])

    # runs along the table and calculates the polynomial according to L.
    for i in range(len(table)):
        l_of_i = 1
        for j in range(len(table)):
            if i != j:
                l_of_i *= (x - x_list[j]) / (x_list[i] - x_list[j])
        p_of_n += l_of_i * y_list[i]

    print(f"The result of the {name} interpolation on {x}  is: p({x}) = {round(p_of_n, 10)}")


def Print_Spline_table(table):
    """
    function that print spline table

    Parameters
    ----------
    table : Points table

    Returns
    -------
    None

    """
    par_name_sp = ["i", "x", "y", "h"]
    table_with_index_sp = [[i, row[0], row[1], round(table[i + 1][0] - row[0], 6) if i != len(table) - 1 else "null"]
                           for
                           row, i in zip(table, range(len(table)))]
    print(tabulate(table_with_index_sp, par_name_sp, 'fancy_grid'))


def Print_table(table):
    """
    function that print table

    Parameters
    ----------
    table : Points table

    Returns
    -------
    None

    """
    par_name = ["i", "x", "y"]
    table_with_index = [[i, row[0], row[1]] for row, i in zip(table, range(len(table)))]
    print(tabulate(table_with_index, par_name, 'fancy_grid'))


def insert_table(table):
    """
    function that insert table- default or by input

    Parameters
    ----------
    table : Points table

    Returns
    -------
    new table

    """
    choice = input("Welcome! Please choose:\n"
                   "1.Get default table\n"
                   "2.Insert self parameters of points\n")
    if choice == "2":
        table.clear()
        n = int(input("Insert amount of points you want to insert"))
        for _ in range(n):
            x = input("x: ")
            y = input("y: ")
            table.append([x, y])

    return table


def MainFunction():

    choice = True
    while choice:
        print("---\n")
        choice = int(input(f"Welcome! We want to find an interpolation with given table x value"
                           "\nChoose Method of solving:\n"
                           "\t1. Linear interpolation.\n"
                           "\t2. Polynomial interpolation.\n"
                           "\t3. Lagrange interpolation.\n"
                           "\t4. Neville interpolation.\n"
                           "\t5. Spline Cubic interpolation\n"
                           "\tEXIT. Click any other char.\n"
                           "---\n"))

        if choice == 1:
            table = [[0, 0], [1, 0.8415], [2, 0.9093], [3, 0.1411], [4, -0.7568], [5, -0.9589]]
            X = 2.5
            linear_interpolation(table, X)
        elif choice == 2:
            table = [[1, 0.8415], [2, 0.9093], [3, 0.1411]]
            X = 2.5
            polynomial_interpolation(table, X)
        elif choice == 3:
            table = [[1, 1], [2, 0], [4, 1.5]]
            X = 3
            lagrange_interpolation(table, X)
        elif choice == 4:
            table = [[1, 0], [1.2, 0.112463], [1.3, 0.167996], [1.4, 0.222709]]
            X = 1.28
            neville_interpolation(table, X)
        elif choice == 5:
            table = [[0, 0], [math.pi / 6, 0.5], [math.pi / 4, 0.7072], [math.pi / 2, 1]]
            X = math.pi / 3
            spline_interpolation(table, X)
        else:
            print("Goodbye :)")
            break

    return


MainFunction()