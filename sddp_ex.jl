using HiGHS
using SDDP

graph = SDDP.LinearGraph(3)

function subproblem_builder(subproblem::Model, node::Int)
    # State variables
    @variable(subproblem, 0 <= volume <= 200, SDDP.State, initial_value = 200)
    # Control variables
    @variables(subproblem, begin
        thermal_generation >= 0
        hydro_generation >= 0
        hydro_spill >= 0
    end)
    # Random variables
    @variable(subproblem, inflow)
    Ω = [0.0, 50.0, 100.0]
    P = [1 / 3, 1 / 3, 1 / 3]
    SDDP.parameterize(subproblem, Ω, P) do ω
        return JuMP.fix(inflow, ω)
    end
    # Transition function and constraints
    @constraints(
        subproblem,
        begin
            volume.out == volume.in - hydro_generation - hydro_spill + inflow
            demand_constraint, hydro_generation + thermal_generation == 150
        end
    )
    # Stage-objective
    if node == 1
        @stageobjective(subproblem, 50 * thermal_generation)
    elseif node == 2
        @stageobjective(subproblem, 100 * thermal_generation)
    else
        @assert node == 3
        @stageobjective(subproblem, 150 * thermal_generation)
    end
    return subproblem
end

model = SDDP.PolicyGraph(
    subproblem_builder,
    graph;
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
)

SDDP.train(model; iteration_limit = 10)

rule = SDDP.DecisionRule(model; node = 1)

solution = SDDP.evaluate(
    rule;
    incoming_state = Dict(:volume => 150.0),
    noise = 50.0,
    controls_to_record = [:hydro_generation, :thermal_generation],
)