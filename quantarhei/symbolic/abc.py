try:
    from sympy import Symbol as _Symbol

    a = _Symbol("a", real=True)
    b = _Symbol("b", real=True)
    c = _Symbol("c", real=True)
    d = _Symbol("d", real=True)
    e = _Symbol("e", real=True)
    f = _Symbol("f", real=True)

    s = _Symbol("s", real=True)
    t = _Symbol("t", real=True)
    T = _Symbol("T", real=True)
    tau = _Symbol("tau", real=True)
    x = _Symbol("x", real=True)
    y = _Symbol("y", real=True)
    t1 = _Symbol("t1", real=True)
    t2 = _Symbol("t2", real=True)
    t3 = _Symbol("t3", real=True)
    tau1 = _Symbol("tau1", real=True)
    tau2 = _Symbol("tau2", real=True)
    tau3 = _Symbol("tau3", real=True)
except ImportError:
    pass
