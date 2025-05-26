import math

def binomial_weight(N, W, p):
    if N < 20:
        return math.comb(N, W) * (p**W) * ((1 - p)**(N - W))
    else:
        lam = N * p
        log_pmf = (-lam) + W * math.log(lam) - math.lgamma(W + 1)
        return math.exp(log_pmf)

# Utility: exact binomial PMF for comparison
def exact_binomial(N, W, p):
    return math.comb(N, W) * (p**W) * ((1 - p)**(N - W))

# Run tests
def run_tests():
    test_cases = [
        # Small N, exact match
        (10, 3, 0.2),
        (10, 0, 0.5),
        (10, 10, 0.9),

        # Medium N
        (50, 10, 0.1),
        (100, 50, 0.5),
        (150, 0, 0.1),

        # Large N (Poisson domain)
        (300, 5, 0.01),        # lam = 3
        (500, 2, 0.004),       # lam = 2
        (1000, 10, 0.01),      # lam = 10
        (1000, 0, 0.01),
        (1000, 20, 0.01),      # check tail accuracy
        (1000, 1000, 0.99),    # extreme case

        # Edge cases
        (100, 0, 0.0),
        (100, 0, 1.0),
        (100, 100, 1.0),
        (100, 100, 0.0)
    ]

    for N, W, p in test_cases:
        result = binomial_weight(N, W, p)
        ref = exact_binomial(N, W, p)


        error = abs(result - ref)
        rel_error = error / (ref + 1e-15)
        print(f"N={N:<5} W={W:<4} p={p:<5}  result={result:.4e}  ref={ref:.4e}  rel_error={rel_error:.2e}")




if __name__ == "__main__":

    run_tests()
