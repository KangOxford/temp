class FirstExact(SinglePatternODESolver):
    r"""
    Solves 1st order exact ordinary differential equations.

    A 1st order differential equation is called exact if it is the total
    differential of a function. That is, the differential equation

    .. math:: P(x, y) \,\partial{}x + Q(x, y) \,\partial{}y = 0

    is exact if there is some function `F(x, y)` such that `P(x, y) =
    \partial{}F/\partial{}x` and `Q(x, y) = \partial{}F/\partial{}y`.  It can
    be shown that a necessary and sufficient condition for a first order ODE
    to be exact is that `\partial{}P/\partial{}y = \partial{}Q/\partial{}x`.
    Then, the solution will be as given below::

        >>> from sympy import Function, Eq, Integral, symbols, pprint
        >>> x, y, t, x0, y0, C1= symbols('x,y,t,x0,y0,C1')
        >>> P, Q, F= map(Function, ['P', 'Q', 'F'])
        >>> pprint(Eq(Eq(F(x, y), Integral(P(t, y), (t, x0, x)) +
        ... Integral(Q(x0, t), (t, y0, y))), C1))
                    x                y
                    /                /
                   |                |
        F(x, y) =  |  P(t, y) dt +  |  Q(x0, t) dt = C1
                   |                |
                  /                /
                  x0               y0

    Where the first partials of `P` and `Q` exist and are continuous in a
    simply connected region.

    A note: SymPy currently has no way to represent inert substitution on an
    expression, so the hint ``1st_exact_Integral`` will return an integral
    with `dy`.  This is supposed to represent the function that you are
    solving for.

    Examples
    ========

    >>> from sympy import Function, dsolve, cos, sin
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> dsolve(cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x),
    ... f(x), hint='1st_exact')
    Eq(x*cos(f(x)) + f(x)**3/3, C1)

    References
    ==========

    - https://en.wikipedia.org/wiki/Exact_differential_equation
    - M. Tenenbaum & H. Pollard, "Ordinary Differential Equations",
      Dover 1963, pp. 73

    # indirect doctest

    """
    hint = "1st_exact"
    has_integral = True
    order = [1]

    def _wilds(self, f, x, order):
        P = Wild('P', exclude=[f(x).diff(x)])
        Q = Wild('Q', exclude=[f(x).diff(x)])
        return P, Q

    def _equation(self, fx, x, order):
        P, Q = self.wilds()
        return P + Q*fx.diff(x)

    def _verify(self, fx) -> bool:
        P, Q = self.wilds()
        x = self.ode_problem.sym
        y = Dummy('y')

        m, n = self.wilds_match()

        m = m.subs(fx, y)
        n = n.subs(fx, y)
        numerator = cancel(m.diff(y) - n.diff(x))

        if numerator.is_zero:
            # Is exact
            return True
        else:
            # The following few conditions try to convert a non-exact
            # differential equation into an exact one.
            # References:
            # 1. Differential equations with applications
            # and historical notes - George E. Simmons
            # 2. https://math.okstate.edu/people/binegar/2233-S99/2233-l12.pdf

            factor_n = cancel(numerator/n)
            factor_m = cancel(-numerator/m)
            if y not in factor_n.free_symbols:
                # If (dP/dy - dQ/dx) / Q = f(x)
                # then exp(integral(f(x))*equation becomes exact
                factor = factor_n
                integration_variable = x
            elif x not in factor_m.free_symbols:
                # If (dP/dy - dQ/dx) / -P = f(y)
                # then exp(integral(f(y))*equation becomes exact
                factor = factor_m
                integration_variable = y
            else:
                # Couldn't convert to exact
                return False

            factor = exp(Integral(factor, integration_variable))
            m *= factor
            n *= factor
            self._wilds_match[P] = m.subs(y, fx)
            self._wilds_match[Q] = n.subs(y, fx)
            return True

    def _get_general_solution(self, *, simplify_flag: bool = True):
        m, n = self.wilds_match()
        fx = self.ode_problem.func
        x = self.ode_problem.sym
        (C1,) = self.ode_problem.get_numbered_constants(num=1)
        y = Dummy('y')

        m = m.subs(fx, y)
        n = n.subs(fx, y)

        gen_sol = Eq(Subs(Integral(m, x)
                          + Integral(n - Integral(m, x).diff(y), y), y, fx), C1)
        return [gen_sol]
