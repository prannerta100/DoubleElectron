def flub(n):
    """Print the Fibonacci series up to n."""
#    a, b = 0, 1
#    while b < n:
#        print(b, end=' ')
#        a, b = b, a + b
    cnt = 0
    for i in range(n):
        for j in range(n):
            for k in range(10):
                cnt += 1
    print()
