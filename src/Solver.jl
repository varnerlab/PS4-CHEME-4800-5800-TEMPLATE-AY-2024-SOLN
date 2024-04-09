function _mysolve(problem::MySimpleProblemModel, solver::MyRungeKuttaMethod)::Tuple{Array{Float64,1}, Array{Float64,2}}
    
    # initialize -
    t0, t1, dt = problem.time_span;
    ic = problem.initial_conditions;

    # get the parameters -
    parameters = problem.parameters;
    p = Array{Any,1}(undef,8)
    p[1] = parameters["α"];
    p[2] = parameters["β"];
    p[3] = parameters["K"];
    p[4] = parameters["K̂"];
    p[5] = parameters["μmax"];
    p[6] = parameters["Y"];
    p[7] = parameters["emax"];
    p[8] = parameters["number_of_enzymes"];
    
    # create the ODEProblem -
    prob = ODEProblem(balances, ic, (t0, t1), p);
    soln = solve(prob, RK4(), reltol=1e-8, abstol=1e-8, saveat=dt);

    # @show typeof(soln)

    # build soln and time arrays -
    T = soln.t
    tmp = soln.u
    number_of_time_steps = length(T)
    number_of_dynamic_states = length(ic);
    X = Array{Float64,2}(undef, number_of_time_steps,  number_of_dynamic_states);
    for i ∈ 1:number_of_time_steps
        soln_vector = tmp[i]
        for j ∈ 1:number_of_dynamic_states
            X[i,j] = soln_vector[j]
        end
    end

    # return - 
    return (T, X)
end

"""
    solve(balances::Function, tspans::Tuple{Float64,Float64,Float64}, initial_conditions::Array{Float64,1}, parameters::Dict{String,Array{Float64,1}}; 
        solver::AbstractIVPSolverType = MyForwardEulerMethod())::Tuple{Array{Float64,1}, Array{Float64,2}}

Solve the system of ODEs defined by the `balances` function, with the given `initial_conditions` and `parameters` over the time span `tspans`. 
The `solver` keyword argument is used to specify the method to solve the system of ODEs. The default solver is the `MyForwardEulerMethod`.

### Arguments
- `balances::Function`: the function that defines the system of ODEs.
- `tspans::Tuple{Float64,Float64,Float64}`: the time span of the simulation.
- `initial_conditions::Array{Float64,1}`: the initial conditions of the system of ODEs.
- `parameters::Dict{String,Array{Float64,1}}`: the parameters of the system of ODEs.

### Optional keyword arguments
- `solver::AbstractIVPSolverType = MyForwardEulerMethod()`: the method to solve the system of ODEs.

### Returns
- `Tuple{Array{Float64,1}, Array{Float64,2}}`: a tuple containing the time array and the solution array.
"""
function mysolve(balances::Function, tspan::Tuple{Float64,Float64,Float64}, initial::Array{Float64,1}, parameters::Dict{String, Any}; 
    solver::AbstractIVPSolverType = MyForwardEulerMethod())::Tuple{Array{Float64,1}, Array{Float64,2}}
    
    # create the problem, object -
    problem = build(MySimpleProblemModel, (
        parameters = parameters,
        initial_conditions = initial,
        time_span = tspan,
        model = balances
    ));

    # solve the problem using the appropriate solver -
    return _mysolve(problem, solver)
end