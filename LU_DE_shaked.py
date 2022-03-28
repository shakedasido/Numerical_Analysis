
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
        elementary[i][i] = 1/mtx[i][i]
        # gets the result out of the multiplication of the elementary matrix with the last result.
        result = mul_mtx_with_mtx(elementary, result)
        mtx = mul_mtx_with_mtx(elementary, mtx)
        # turn to zero all the elements below the pivot
        for j in range(i+1, len(mtx)):
            elementary = create_mtx_i(len(mtx[0]), len(mtx))
            elementary[j][i] = -(mtx[j][i])
            mtx = mul_mtx_with_mtx(elementary, mtx)
            result = mul_mtx_with_mtx(elementary, result)
    # turn to zero all the elements above the pivot
    for i in range(len(mtx[0])-1, 0, -1):
        for j in range(i-1, -1, -1):
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
    print_mul(mtx1, mtx2, new_matrix)
    print("\n")
    return new_matrix


def print_mul(mtx1, mtx2, res):
    """
        gets two matrix's and print the multiplication.

        Parameters
        ----------
        mtx1 : nested list (nXn)
        mtx2 : nested list (nXn)
        res: nested list (nXn)

        Returns
        -------
        None
    """
    for i in range(0, len(mtx1)):
        if i != int(len(mtx1) / 2):
            print("[" + '  '.join(map(str, mtx1[i])) + "]" + "    [" + '  '.join(
                map(str, mtx2[i])) + "]" + "  [" + '  '.join(map(str, res[i])) + "]")
        else:
            print("[" + '  '.join(map(str, mtx1[i])) + "]" + " X [" + '  '.join(
                map(str, mtx2[i])) + "]" + "= [" + '  '.join(map(str, res[i])) + "]")


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
                # mul the matrix with the elementary matrix and print it
                print("\nAfter the pivoting process\n")
                mtx = mul_mtx_with_mtx(elementary_switched, mtx)
                vector_b = mul_mtx_with_mtx(elementary_switched, vector_b)
                pivot = abs(mtx[i][i])

    print("The matrix after changing the rows:")
    print_matrix(mtx)
    print("\n")
    print_matrix(vector_b)
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


# MAIN:

mtxA = [[1, 0, -1], [-0.5, 1, -0.25], [1, -0.5, 1]]
vector_b = [[0.2], [-1.425], [2]]
print("\nMatrix A: \n")
print_matrix(mtxA)
print("\nVector b: \n")
print_matrix(vector_b)
# let's switch rows in order to get the bigger pivot
mtxA, vector_b = switch_rows(mtxA, vector_b)
print("\ncreate matrix U: \n")
u_mtx, l_mtx = u_and_l_matrix(mtxA, vector_b)
# U matrix
print("\nMatrix U: \n")
print_matrix(u_mtx)
# L matrix
print("\nMatrix L: \n")
print_matrix(l_mtx)

# solving the matrix and get x, the result.

print("\nU(-)*L(-1)*VECTOR_B = X: \n")
print("\nCreate U^(-1):\n")
U_inverse = get_matrix_inverse(u_mtx, vector_b)
print("\nU^(-1):\n")
print_matrix(U_inverse)
print("\nCreate L^(-1):\n")
L_inverse = get_matrix_inverse(l_mtx, vector_b)
print("\nL^(-1):\n")
print_matrix(L_inverse)
print("\nCreate U^(-1)*l^(-1)* = x: \n")
x = mul_mtx_with_mtx(U_inverse, mul_mtx_with_mtx(L_inverse, vector_b))

print("\nMatrix solve\nx: ")
print_matrix(x)
