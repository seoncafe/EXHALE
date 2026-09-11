"""Small algebraic checks for the solver-approach review.

Run with Python 3. No EXHALE state is loaded or evolved. These examples test
general mathematical implications, not the convergence of a planet model.
The program writes no files and requires only the standard library.
"""

from math import exp, pi, sin, sqrt


def fixed_grid_residual():
    print("1. Centered Poisson equation: algebraic residual vs. solution error")
    for n in (31, 63, 127):
        h = 1.0 / (n + 1)
        eigenvalue = 4.0 * sin(pi * h / 2.0) ** 2 / h**2
        exact = [sin(pi * j * h) for j in range(n + 2)]
        exact[0] = exact[-1] = 0.0
        discrete = [pi**2 / eigenvalue * value for value in exact]
        residual = max(
            abs((2 * discrete[j] - discrete[j - 1] - discrete[j + 1]) / h**2
                - pi**2 * exact[j])
            for j in range(1, n + 1)
        ) / pi**2
        error = max(abs(a - b) for a, b in zip(exact, discrete))
        assert residual < 1e-9 and error > 1e-5
        print(f"   N={n:3d}: normalized residual={residual:.6e}, error={error:.6e}")


def similarity_and_gmres():
    print("2. Similar matrices [[1,a],[0,2]], b=(0,1): one GMRES step")
    for a in (100.0, 1.0):
        alpha = 2.0 / (a * a + 4.0)
        residual = sqrt((alpha * a) ** 2 + (1.0 - 2 * alpha) ** 2)
        print(f"   a={a:g}: eigenvalues=(1,2), relative residual={residual:.9f}")
    assert sqrt(1 - 4 / 10004) > 0.99
    assert sqrt(1 - 4 / 5) < 0.45


def compression_and_coupling():
    print("3. J=[[1,1],[1,0]]: species compression is zero")
    condition = (3 + sqrt(5)) / 2
    print(f"   full determinant=-1, full 2-norm condition={condition:.9f}, Schur=-1")
    print("4. J=[[1,2],[2,1]]: exact scalar block solves need not converge jointly")
    for omega in (1.0, 0.5, 0.125):
        gain = 1 - omega + 4 * omega
        assert gain > 1
        print(f"   omega={omega:g}: composition iteration gain={gain:g}")


def split_steady_state():
    print("5. Exact split flows x'=1-x, then x'=-x; full steady state=0.5")
    for dt in (0.1, 0.05, 0.025):
        x = 1 / (exp(dt) + 1)
        split_x = exp(-dt) * (1 + (x - 1) * exp(-dt))
        assert abs(split_x - x) < 1e-15
        residual = abs(1 - 2 * x)
        assert residual > 0.01
        print(f"   dt={dt:g}: split fixed point={x:.9f}, unsplit residual={residual:.9f}")


def residual_as_its_own_scale():
    print("6. Residual used as its own reference scale")
    for residual in (1e-4, 1e-8, 1e-12):
        scaled = abs(residual) / max(abs(residual), 1e-30)
        assert scaled == 1.0
        print(f"   residual={residual:.0e}: normalized residual={scaled:g}")


if __name__ == "__main__":
    fixed_grid_residual()
    similarity_and_gmres()
    compression_and_coupling()
    split_steady_state()
    residual_as_its_own_scale()
    print("All illustrative assertions passed; no EXHALE regression was run.")
