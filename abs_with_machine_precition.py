
epsilon = 1.0
flag = True
while flag:
    final_epsilon = epsilon
    epsilon = epsilon / 2.0
    if 1.0 + epsilon == 1.0:
        flag = False

print(abs(3.0 * (4.0 / 3.0 - 1) - 1) - final_epsilon)
