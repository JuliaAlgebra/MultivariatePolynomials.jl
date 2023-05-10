using RowEchelon
n = 2
m = 6
M =
    [
        1.0000 1.5868 2.2477 2.7603 3.6690 5.2387
        1.5868 2.7603 3.6690 5.1073 6.5115 8.8245
        2.2477 3.6690 5.2387 6.5115 8.8245 12.7072
        2.7603 5.1073 6.5115 9.8013 12.1965 15.9960
        3.6690 6.5115 8.8245 12.1965 15.9960 22.1084
        5.2387 8.8245 12.7072 15.9960 22.1084 32.1036
    ] + 1.44e-5 * eye(m)

#d, U = eigs(M, nev=3)
#d, U = eigs(M)
#function getU(U)
#    V = [U[:,1]*sqrt(d[1]) U[:,2]*sqrt(d[2]) U[:,3]*sqrt(d[3])]
#    Matrix{Int}(round.(rref(V')'))
#end
r = rank(M, 1e-4)
V = chol(M)[1:r, :]
U = rref(V)'
Ns = [
    [
        0 1 0
        -2 3 0
        -4 2 2
    ],
    [
        0 0 1
        -4 2 2
        -6 0 5
    ],
]
#N = 0.6906Ns[1] + 0.3091Ns[2]
for i in 1:10
    λ = rand()
    N = λ * Ns[1] + (1 - λ) * Ns[2]
    Z = schurfact(N)[:Z]
    x = [Vector{Float64}(n) for j in 1:r]
    for j in 1:r
        qj = Z[:, j]
        for i in 1:n
            x[j][i] = dot(qj, Ns[i] * qj)
        end
    end
    @show x
end
