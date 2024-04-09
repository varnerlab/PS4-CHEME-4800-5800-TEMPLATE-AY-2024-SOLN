function _compute_enzyme_synthesis_rates(x,p)

    # order of the parameters 
    # p[1] = α
    # p[2] = β
    # p[3] = K
    # p[4] = K̂
    # p[5] = μmax
    # p[6] = Y
    # p[7] = emax
    # p[8] = 2; # number of enzymes


    # get parameters -
    α = p[1]; # enzyme synthesis rate constants
    K̂ = p[4]; # saturation constants for enzyme synthesis
    number_of_enzymes = p[8];

    # initialize the synthesis rates vector -
    rE = zeros(number_of_enzymes);

    # alias the state variables -
    S₁ = x[1];
    S₂ = x[2];

    # compute the synthesis rates -
    rE[1] = α[1]*S₁/(K̂[1] + S₁);
    rE[2] = α[2]*S₂/(K̂[2] + S₂);

    # return the synthesis rates -
    return rE;
end

function _compute_cybernetic_variables(rG, p)

    # get parameters -
    number_of_enzymes = p[8];
    u = zeros(number_of_enzymes);
    v = zeros(number_of_enzymes);

    # compute the u-variable -
    u[1] = rG[1]/(sum(rG));
    u[2] = rG[2]/(sum(rG));

    # compute the v-variable -
    v[1] = rG[1]/maximum(rG);
    v[2] = rG[2]/maximum(rG);

    # return the cybernetic variables -
    return (u,v);
end

function _compute_growth_rates(x,p)

    # order of the parameters 
    # p[1] = α
    # p[2] = β
    # p[3] = K
    # p[4] = K̂
    # p[5] = μmax
    # p[6] = Y
    # p[7] = emax
    # p[8] = 2; # number of enzymes

    # get parameters -
    μmax = p[5]; # max growth rate
    emax = p[7]; # max enzyme concentration
    K = p[3];    # saturation constants for growth
    number_of_enzymes = p[8];
    
    # alias the state variables -
    S₁ = x[1];
    S₂ = x[2];
    e₁ = x[3];
    e₂ = x[4];

    # compute the growth rates -
    rG = zeros(number_of_enzymes);

    # compute the growth rates -
    rG[1] = μmax[1]*(e₁/emax[1])*(S₁/(K[1] + S₁));
    rG[2] = μmax[2]*(e₂/emax[2])*(S₂/(K[2] + S₂));

    # return the growth rates -
    return rG;
end


# --- PUBLIC METHODS BELOW HERE ----------------------------------------------------------------------------------- #
"""
    balances(dx, x, p, t)

Computes the right-hand sides of the balances for the Kompala model, Biotechnol. Bioeng. 1986, 1044-1055. 

### Arguments
- `dx::Array{Float64,1}`: the right-hand-side (rhs) of the balances
- `x::Array{Float64,1}`: the state variables
- `p::Array{Float64,1}`: the parameters
- `t::Float64`: the time variable
"""
function balances(dx, x, p, t)

    # order of the parameters 
    # p[1] = α
    # p[2] = β
    # p[3] = K
    # p[4] = K̂
    # p[5] = μmax
    # p[6] = Y
    # p[7] = emax
    # p[8] = 2; # number of enzymes

    # initialize -
    Y = p[6];    # yield coefficients
    β = p[2];    # decay rate constant
    e₁, e₂, C = x[3], x[4], x[5];  # alias the state variables -

    # compute the kinetics
    rE = _compute_enzyme_synthesis_rates(x,p);
    rG = _compute_growth_rates(x,p);
    (u,v) = _compute_cybernetic_variables(rG, p);  # compute the control (cybernetic variables) -

    # compute the dilution term -
    μ = sum(rG.*v);

    # setup the right-hand-side (rhs) -
    dx[1] = -(1/Y[1])*rG[1]*v[1]*C;
    dx[2] = -(1/Y[2])*rG[2]*v[2]*C;
    dx[3] = rE[1]*u[1] - (μ + β[1])*e₁;
    dx[4] = rE[2]*u[2] - (μ + β[2])*e₂;
    dx[5] = μ*C;
end
# --- PUBLIC METHODS ABOVE HERE ----------------------------------------------------------------------------------- #