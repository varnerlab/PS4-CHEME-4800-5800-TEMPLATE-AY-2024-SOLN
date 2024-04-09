# setup paths -
const _ROOT = @__DIR__;
const _PATH_TO_SRC = joinpath(_ROOT, "src");
const _PATH_TO_TEST_DATA = joinpath(_ROOT, "testdata");
const _PATH_TO_MY_DATA = joinpath(_ROOT, "mydata");

# update packages -
using Pkg;
Pkg.activate("."); Pkg.resolve(); Pkg.instantiate(); Pkg.update();

# load external packages -
using DifferentialEquations;
using JLD2
using FileIO
using Test


# include files -
include(joinpath(_PATH_TO_SRC, "Types.jl"));
include(joinpath(_PATH_TO_SRC, "Factory.jl"));
include(joinpath(_PATH_TO_SRC, "Kompala.jl"));
include(joinpath(_PATH_TO_SRC, "Solver.jl"));


