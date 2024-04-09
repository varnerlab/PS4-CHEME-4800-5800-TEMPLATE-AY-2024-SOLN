# load the price dataset full dataset, remove firms with missing data, 
# compute the growth rate matrix, and then compute the singular value decomposition of the growth rate matrix.
include("Include.jl") 

# ----------------------------------------------------------------------------------
# for more information on tests, see: https://docs.julialang.org/en/v1/stdlib/Test/
# ----------------------------------------------------------------------------------

# Testset - let's write some tests for all the facts that we know about the SVD -
@testset verbose = true "Test the Kompala Simulation Refactoring" begin

    # load the test data, and your simulation data
    mydata = load(joinpath(_PATH_TO_MY_DATA, "mysimulationdata.jld2"));
    testdata = load(joinpath(_PATH_TO_TEST_DATA, "testdata.jld2"));

    # get the components from the test data -
    T1 = mydata["T"];
    X1 = mydata["X"];
    T2 = testdata["T"];
    X2 = testdata["X"];

    # test the time arrays -
    @test length(T1) == length(T2);
    @test round.(T1, digits=4) ≈ round.(T2, digits=4);

    # test the state arrays -
    @test size(X1) == size(X2);
    @test round.(X1, digits=4) ≈ round.(X2, digits=4);
end