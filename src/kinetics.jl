# Equilibria of reaction schemes, for the initial conditions of kinetic models in model_zoo/.

using LinearAlgebra

"""
    equilibrium(A, total) → Vector
Equilibrium of the closed linear system dx/dt = A x whose states add up to `total`. The rows
of A x = 0 are linearly dependent, so the last state's row is replaced by the conservation
sum(x) = total and the system is solved.
"""
function equilibrium(A, total)
    M = [A[1:(end - 1), :]; ones(eltype(A), 1, size(A, 2))]
    b = [zeros(eltype(A), size(A, 1) - 1); total]
    return M \ b
end
