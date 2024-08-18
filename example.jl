using GLPK
using HiGHS
using StochasticPrograms

@stochastic_model simple_model begin
    @stage 1 begin
        @decision(simple_model, x >= 0)
        @objective(simple_model, Max, -30*x)
    end
    @stage 2 begin
        @uncertain d
        @recourse(simple_model, 0 <= s)
        @recourse(simple_model, 0 <= i)
        @recourse(simple_model, 0 <= l)
        @objective(simple_model, Max, 60*s - 10*i - 5*l)
        @constraint(simple_model, d == s + l)
        @constraint(simple_model, i == x - s)
    end
end

e1 = @scenario d = 45 probability = 0.7
e2 = @scenario d = 40 probability = 0.2
e3 = @scenario d = 50 probability = 0.1

sp = instantiate(simple_model, [e1, e2, e3], optimizer = GLPK.Optimizer)

optimize!(sp)

obj = objective_value(sp)