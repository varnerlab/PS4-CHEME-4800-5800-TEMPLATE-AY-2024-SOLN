# include -
include("Include.jl");

# setup parameters for the model -
α = [1e-3,1e-3];
β = [0.05,0.05];
K = [0.01,0.2];
K̂ = [0.01,0.2];
μmax = [1.08,0.82];
Y = [0.52,0.5];
emax = [
    α[1]/(μmax[1] + β[1]),
    α[2]/(μmax[2] + β[2])
];

# setup the initial conditions -
xₒ = [
    0.5;            # 1 S₁
    2.5;            # 2 S₂
    0.90*emax[1];   # 3 E₁
    0.18*emax[2];   # 4 E₂
    4e-3;           # 5 C
];

parameters = Dict{String,Any}(
    "α" => α,
    "β" => β,
    "K" => K,
    "K̂" => K̂,
    "μmax" => μmax,
    "Y" => Y,
    "emax" => emax,
    "number_of_enzymes" => 2
);

# call the solver -
(T,X) = mysolve(balances, (0.0, 10.0, 0.1), xₒ, parameters, solver = MyRungeKuttaMethod());

# save the results in a dictionary -
save(joinpath(_PATH_TO_MY_DATA, "mysimulationdata.jld2"), Dict("T" => T, "X" => X));