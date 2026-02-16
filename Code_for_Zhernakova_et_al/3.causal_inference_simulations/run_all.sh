#!/bin/bash
# Continue execution even if a command fails
set +e

# ----------------------------
# Run linear and quadratic simulations for different beta values
# ----------------------------
for beta in -1.2 -0.2 0 0.2 1.2; do
    echo "Running simulations for beta = $beta"
    Rscript simulations_linear.R "$beta" || echo "Failed: simulations_linear.R $beta"
    Rscript simulations_quadratic.R "$beta" || echo "Failed: simulations_quadratic.R $beta"
done

# ----------------------------
# Run missing data simulations
# ----------------------------
for perc_miss in 0.1 0.2 0.3; do
    for miss_type in MNAR MAR; do
        for rem in TRUE FALSE; do
            echo "Running missing data simulation: perc_miss=$perc_miss, type=$miss_type, remove=$rem"
            Rscript Final_parallelized_missing.R "$perc_miss" "$miss_type" "$rem" || \
                echo "Failed: Final_parallelized_missing.R $perc_miss $miss_type $rem"
        done
    done
done

# ----------------------------
# Run non-linear simulations
# ----------------------------
for rel in cos cubic exp linear log quadratic; do
    echo "Running non-linear simulations for relation: $rel"
    Rscript simulations_nonLinear.R "$rel" || echo "Failed: simulations_nonLinear.R $rel"
    Rscript simulations_nonLinear_with_t.R "$rel" || echo "Failed: simulations_nonLinear_with_t.R $rel"
    Rscript Final_linear_on_non_linear.R "$rel" || echo "Failed: Final_linear_on_non_linear.R $rel"
    Rscript Final_linear_on_non_linear_with_t.R "$rel" || echo "Failed: Final_linear_on_non_linear_with_t.R $rel"
done

# ----------------------------
# Run covariance simulations
# ----------------------------
echo "Running covariance simulations"
Rscript simulations_cov_testing.R || echo "Failed: simulations_cov_testing.R"
Rscript simulations_cov_without_testing.R || echo "Failed: simulations_cov_without_testing.R"
Rscript simulations_cov_larger_testing.R || echo "Failed: simulations_cov_larger_testing.R"
Rscript simulations_cov_larger_without_testing.R || echo "Failed: simulations_cov_larger_without_testing.R"

