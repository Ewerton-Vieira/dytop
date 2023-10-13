def iterate(f, x, n=1):
    i = 1
    y = f(x)
    while i<n:
        y = f(y)
        i += 1
    return y