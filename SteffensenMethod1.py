import matplotlib.pyplot as plt
import numpy as np


# Definition of the Steffensen method:
def SteffensenMethod(p0, f):
    epsilon = 0.000001
    maxNumOfIter = 10000

    for iter in range(1, maxNumOfIter):

        p1 = p0 + f(p0)
        p2 = p1 + f(p1)
        p = p2 - (pow((p2 - p1), 2) / (p2 - (2 * p1) + p0))

        # print the iteration's values:
        print("#", iter, " x= ", p0, "\t\tF(x)= ", f(p0))

        # if the current x, minus the last p0 is lass then the error evaluation, then we have the result:
        if abs(p - p0) < epsilon:
            print(f"\nThe result of x in the {iter} iteration is: {round(p, 4)}\n")
            return

        p0 = p


def main():
    # Steffensen, by given a starting point:
    print("\nSteffensen Method:\n")
    print("---------------------------------\n")

    choice = int(input("\nPlease choose a function: \n"
                       "1. x ^3 - x - 1\n"
                       "2. 2 * x ^3 - 2 * x - 5\n"
                       "3. x ^3 + 2 * x ^2 + x - 1\n"))

    x = 0
    if choice == 1:
        while x != 2 and x != 1.5 and x != 3: x = float(input("Please choose a start point:\n "
                                                              "\t1.5. x_0 = 1.5\n"
                                                              "\t2. x_0 = 2\n"
                                                              "\t3. x_0 = 3\n"
                                                              "---\n"))
        plt.title("x ^3 - x - 1")
        print(f"\nSteffensen with function: f = x ^3 - x - 1, and Start point of: x = {x}\n")

        # Definition of f(x):
        def f(a):
            return a ** 3 - a - 1

        func = f

    elif choice == 2:
        while x != 1 and x != 1.5 and x != 2 and x != 3: x = float(input("Please choose a start point:\n "
                                                                         "\t1. x_0 = 1\n"
                                                                         "\t1.5. x_0 = 1.5\n"
                                                                         "\t2. x_0 = 2\n"
                                                                         "\t3. x_0 = 3\n"
                                                                         "---\n"))
        plt.title("2 * x ^3 - 2 * x - 5")
        print(f"\nSteffensen with function: f = 2 * x ^3 - 2 * x - 5, and Start point of: x = {x}\n")

        # Definition of f(x):
        def f(a):
            return 2 * a ** 3 - 2 * a - 5

        func = f


    elif choice == 3:
        while x != 1 and x != 0.5 and x != 2 and x != 3: x = float(input("Please choose a start point:\n "
                                                                         "\t0.5. x_0 = 0.5\n"
                                                                         "\t2. x_0 = 2\n"
                                                                         "\t3. x_0 = 3\n"
                                                                         "---\n"))
        plt.title("x ^3 + 2 * x ^2 + x - 1")
        print(f"\nSteffensen with function: f = x ^3 + 2 * x ^2 + x - 1, and Start point of: x = {x}\n")

        # Definition of f(x):
        def f(a):
            return a ** 3 + 2 * a ** 2 + a - 1

        func = f

    else:
        print("Goodbye :)")
        return

    SteffensenMethod(x, func)
    y = np.linspace(-2, 2)  # printing range
    plt.plot(y, func(y))
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.grid()
    plt.show()


main()
