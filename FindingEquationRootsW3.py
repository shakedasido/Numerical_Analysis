from math import e
import sympy as sp
import math

EPSILON4 = 0.00001
EPSILON5 = 0.000001
ITER = 10
x = sp.symbols('x')
FUNCTION = x ** 4 + x ** 3 - 3 * x ** 2
checkRange = range(-3, 2)
TOP_ITER = 100
ROUND = 5


def BisectionMethod(XPlacingInFunction, a, b, iteration_count, epsilon):
    """
        Function that find a root by using the Bisection Method by given two guesses and
        By crossing each range into two halves and seeing where the point is by the function sign. Note: a<=b

        Attributes
        ----------
        XPlacingInFunction: func
            The function which we calculate it's roots by using the Bisection Method
        a: float
            name of the start of the range.
        b: float
            name of the end of the range.
        iteration_count: int
            Limit number of iteration
        epsilon: float
            The machine epsilon of the solution

        Return
        ----------
        List unique roots of the equation.


    """
    equationRoots = []
    average = (a + b) / 2

    exp = pow(10, -10)
    if b - a == 0:
        c = 100
    c = ((-1) * math.log(exp / (b - a), e)) / math.log(2, e)

    if c < iteration_count:
        print("The Method cannot converge.")
        return equationRoots

    if (abs(b - a)) < epsilon:
        print("\n** Current range **\n\t", round(a, 1), "to", round(a + 0.1, 1))
        print("The root is: ", round(average, ROUND))
        equationRoots.append(round(average, ROUND))
        return equationRoots

    iteration_count += 1
    if XPlacingInFunction(a) * XPlacingInFunction(average) > 0:
        equationRoots += BisectionMethod(XPlacingInFunction, average, b, iteration_count, epsilon)
        return equationRoots
    else:
        equationRoots += BisectionMethod(XPlacingInFunction, a, average, iteration_count, epsilon)
        return equationRoots


def Bisection(xPlacingInFunction, range_check, epsilon):
    """
        Function that find a root by using the Bisection Method by given two guesses and
        By crossing each range into two halves and seeing where the point is by the function sign. Note: a<=b

        Attributes
        ----------
        xPlacingInFunction: func
            The function which we calculate it's roots by using the Bisection Method
        range_check : List
            guesses of a points that are getting close to the desired root.
        epsilon: float
            The machine epsilon of the solution

    """
    iteration_count = 0
    equationRoots = []

    # first number in the range (int)
    for i in range_check:
        # add 10 times 0.1 to the number, till it reaches the next number in the range.
        for j in range(10):
            a = round((j * 0.1) + i, 2)
            b = round(a + 0.1, 2)

            if xPlacingInFunction(b) == 0:
                print("\n** Current range **\n\t", round(a, 1), "to", round(b, 1))
                print("root in ", b)
                equationRoots.append(b)
            if (xPlacingInFunction(a) * xPlacingInFunction(b)) < 0:
                equationRoots += BisectionMethod(xPlacingInFunction, a, b, iteration_count, epsilon)

    return equationRoots


def NewtonRaphsonMethod(func, point_guess, Attempts, epsilon):
    """
        Function that find a root by using the Newton Raphson Method by given a guess of a point
        that close to the desired root.

        Attributes
        ----------
        func : Function
            The function which we calculate it's roots by using the Secant Method
        point_guess : int
            The guess of a point, close to the desired root.
        Attempts : int
            number of iteration
        epsilon: float
            The machine epsilon of the solution

    """

    if func.subs(x, point_guess) == 0:
        return 0
    for i in range(Attempts):
        print(f"** {i + 1} iteration **\n x: {point_guess}\n f({point_guess}) = {func} ="
              f" {round(func.subs(x, point_guess), 2)} \n"
              f" f({point_guess}) = {sp.diff(func, x)} = {round(sp.diff(func, x).subs(x, point_guess), 2)}\n")

        if sp.diff(func, x).subs(x, point_guess) == 0.0:
            continue

        next_x = (point_guess - func.subs(x, point_guess) / sp.diff(func, x).subs(x, point_guess))

        approx = abs(next_x - point_guess)
        if approx < epsilon:
            print("Root Solution => X = ", round(next_x, 8))
            return next_x
        point_guess = next_x
    print("No solution to be found!")
    return None


def NewtonRaphson(func, range_check, Attempts, epsilon):
    """
            Function that find a root by using the Newton Raphson Method by given a guess of a point
            that close to the desired root.

            Attributes
            ----------
            func : Function
                The function which we calculate it's roots by using the Newton Raphson Method
            range_check : List
                guesses of a points that are getting close to the desired root.
            Attempts : int
                number of total iterations.
            epsilon: float
                The machine epsilon of the solution


    """

    equationRoots = []
    # first number in the range (int)
    for i in range_check:
        guess = i
        # if its the last element in the list
        if i == range_check[-1]:
            break
        for _ in range(10):
            guess += round(0.1, 1)

            root = NewtonRaphsonMethod(func, guess, Attempts, epsilon)

            if root is not None and round(root, ROUND) not in equationRoots:
                equationRoots.append(round(root, ROUND))
            else:
                print("Root is already exist in the equation roots list.")

    return equationRoots


def SecantMethod(XPlacingInFunction, a, b, iteration_count, epsilon):
    """
        Recursive function that find a root by using the Secant Method by given two guesses. Note: a<=b

        Attributes
        ----------
        XPlacingInFunction: func
            The function which we calculate it's roots by using the Secant Method
        a: float
            name of the start of the range.
        b: float
            name of the end of the range.
        iteration_count: int
            Limit number of iteration
        epsilon: float
            The machine epsilon of the solution

        Return
        ----------
        The method itself (its a recursive func) with new data / None

    """

    # reach into the limit possible iteration
    if TOP_ITER < iteration_count:
        return

    # recreation function- if the sub between a to b is smaller then the defined epsilon, then return this last root
    if epsilon > abs(b - a):
        print("The root is: ", round(b, ROUND))
        return round(b, ROUND)

    # a*f(b) - b*f(a) / f(b) - f(a)
    next_b = (a * XPlacingInFunction(b) - b * XPlacingInFunction(a)) / (
            XPlacingInFunction(b) - XPlacingInFunction(a))

    return SecantMethod(XPlacingInFunction, b, next_b, iteration_count + 1, epsilon)


def Secant(func, range_check, epsilon):
    """
        Function that find a root by using the Secant Method by given a range and by given two guesses, it iterates
        and decrease the range until it's close to the desired root. Note: a<=b

        Attributes
        ----------
        func : Function
            The function which we calculate it's roots by using the Secant Method
        range_check : List
            guesses of a points that are getting close to the desired root.
        epsilon: float
            The machine epsilon of the solution

    """
    iteration_count = 0
    equationRoots = []

    # a number in the range (int)
    for i in range_check:
        b = i
        # if its the b element in the list
        if i == range_check[-1]:
            break

        # add 10 times 0.1 to the number, till it reaches the next number in the range.
        for _ in range(10):
            a = round(b, 1)
            b += round(0.1, 1)

            print("\n** Current range **\n\t", round(a, 1), "to", round(b, 1))

            root = SecantMethod(func, a, b, iteration_count, epsilon)
            if root is not None and round(root, ROUND) not in equationRoots:
                equationRoots.append(round(root, ROUND))
            else:
                print("Root is already exist in the equation roots list.")

    return equationRoots


def Main():
    """
        main function. operates all the other method function.
    """
    f = FUNCTION
    equationRoots = []

    def XPlacingInFunction(a):
        return f.subs(x, a)

    f_derivatives = lambda a: sp.diff(f, x).subs(x, a)

    choice = True
    while choice:
        print("---\n")
        print(f"Welcome! We want to find Roots of the equation f(x) = {f} "
              "\nChoose Method of solving:\n"
              "\t1. Bisection Method.\n"
              "\t2. Newton Raphson.\n"
              "\t3. Secant Method.\n"
              "\tEXIT. Click any other char.\n"
              "---\n")
        choice = input("")

        if choice == "1":

            print(f"Multiplicity Roots check:")

            print(f"\nChecks the intersection of {f} with the X-axis:")
            equationRoots += Bisection(XPlacingInFunction, checkRange, EPSILON4)

            print(f"\nChecks for lunch point of {f} with the X-axis, by the function derivative:")
            root = Bisection(f_derivatives, checkRange, EPSILON4)

            # Crosses the x-axis
            if XPlacingInFunction(root) == 0:
                equationRoots += root
                print(equationRoots)
            else:
                print(" Not the root of the equation ")

            if equationRoots != 0:
                print("\nThere are ", len(equationRoots), "equationRoots: ", equationRoots)
            equationRoots.clear()

        elif choice == "2":
            equationRoots += NewtonRaphson(f, checkRange, ITER, EPSILON5)
            if equationRoots != 0:
                print("\nThere are ", len(equationRoots), "equationRoots: ", equationRoots)

            equationRoots.clear()

        elif choice == "3":
            equationRoots += Secant(XPlacingInFunction, checkRange, EPSILON5)
            if equationRoots != 0:
                print("\nThere are", len(equationRoots), "equationRoots:", equationRoots)

            equationRoots.clear()
        else:
            print("Goodbye :)")
            break


Main()
