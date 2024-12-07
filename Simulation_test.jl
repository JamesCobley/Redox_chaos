# Sumulation Test
# Author: [James Cobley]
# Description: Models a redox proteoform system, computes the Lyapunov exponent,
# generates plots, and saves data files.

using Plots, Random, CSV, DataFrames, Statistics

# Set the working directory (default to the current directory)
output_dir = pwd()

# Initialize system state
function initialize_state(r::Int, initial_proteoform::String)
    proteoforms = [bitstring(i) |> s -> lpad(s, r, '0')[end-r+1:end] for i in 0:(2^r-1)]
    state = Dict(pf => 0.0 for pf in proteoforms)
    state[initial_proteoform] = 1.0
    return state, proteoforms
end

# Find allowed transitions
function find_allowed_transitions(pf::String, proteoforms::Vector{String}, r::Int)
    allowed = [pf]
    struct_vec = collect(pf)
    current_k = count(c -> c == '1', struct_vec)

    for i in 1:r
        new_structure = copy(struct_vec)
        new_structure[i] = (new_structure[i] == '0' ? '1' : '0')
        new_pf = join(new_structure)
        new_k = count(c -> c == '1', new_structure)
        if abs(new_k - current_k) == 1
            push!(allowed, new_pf)
        end
    end
    return allowed
end

# Evolve the system probabilistically
function evolve_state(state::Dict{String, Float64}, proteoforms::Vector{String}, r::Int, p_jump::Float64)
    new_state = Dict(pf => 0.0 for pf in proteoforms)
    for pf in proteoforms
        current_prob = state[pf]
        allowed_transitions = find_allowed_transitions(pf, proteoforms, r)
        num_allowed = length(allowed_transitions) - 1

        for new_pf in allowed_transitions
            if new_pf != pf && num_allowed > 0
                new_state[new_pf] += current_prob * p_jump / num_allowed
            else
                new_state[pf] += current_prob * (1 - p_jump)
            end
        end
    end
    return new_state
end

# Calculate entropy and mean oxidation state
function calc_metrics(state::Dict{String, Float64}, proteoforms::Vector{String}, r::Int)
    total_molecules = sum(values(state))
    probabilities = [state[pf] / total_molecules for pf in proteoforms]
    entropy = -sum(p * log(p + 1e-10) for p in probabilities)
    mean_k = sum(count('1', pf) * state[pf] for pf in proteoforms) / total_molecules
    return entropy, mean_k
end

# Corrected Lyapunov exponent calculation
function compute_lyapunov_exponent(
    r::Int, initial_proteoform::String, steps::Int, p_jump::Float64, epsilon::Float64
)
    original_state, proteoforms = initialize_state(r, initial_proteoform)
    perturbed_state = deepcopy(original_state)
    perturbed_state[proteoforms[2]] += epsilon
    perturbed_state[proteoforms[1]] -= epsilon

    distances = []
    
    for _ in 1:steps
        dist = sqrt(sum((original_state[pf] - perturbed_state[pf])^2 for pf in proteoforms))
        push!(distances, dist)
        
        original_state = evolve_state(original_state, proteoforms, r, p_jump)
        perturbed_state = evolve_state(perturbed_state, proteoforms, r, p_jump)
        
        perturbed_total = sum(values(perturbed_state))
        for pf in proteoforms
            perturbed_state[pf] /= perturbed_total
        end
    end
    
    lyapunov_exponent = mean(log.(distances .+ 1e-10))
    return lyapunov_exponent
end

# Save Data
function save_history(history::Vector{Dict{String, Float64}}, proteoforms::Vector{String})
    df = DataFrame(TimeStep=Int[], Proteoform=String[], Count=Float64[])
    for (t, state) in enumerate(history)
        for pf in proteoforms
            push!(df, (t, pf, state[pf]))
        end
    end
    CSV.write(joinpath(output_dir, "simulation_history.csv"), df)
end

function save_metrics(steps, entropies, mean_oxidation_states)
    metrics_df = DataFrame(
        TimeStep=1:steps,
        Entropy=entropies,
        MeanOxidationState=mean_oxidation_states
    )
    CSV.write(joinpath(output_dir, "metrics.csv"), metrics_df)
end

# Generate Plots
function plot_poincare(history::Vector{Dict{String, Float64}}, proteoforms::Vector{String}, period::Int)
    points_x = []
    points_y = []
    for t in 1:period:length(history)-period
        x = sum(count('1', pf) * history[t][pf] for pf in proteoforms)
        y = sum(count('1', pf) * history[t+period][pf] for pf in proteoforms)
        push!(points_x, x)
        push!(points_y, y)
    end
    scatter(points_x, points_y, title="Poincaré Diagram", xlabel="k(t)", ylabel="k(t + T)")
    savefig(joinpath(output_dir, "poincare_diagram.png"))
end

function plot_bifurcation(r::Int, initial_proteoform::String, steps::Int, p_min::Float64, p_max::Float64, p_steps::Int)
    p_values = range(p_min, p_max, length=p_steps)
    bifurcation_x = []
    bifurcation_y = []

    for p_jump in p_values
        _, proteoforms, _, mean_k_values = simulate(r, initial_proteoform, steps, p_jump)
        for k in mean_k_values[end-100:end]
            push!(bifurcation_x, p_jump)
            push!(bifurcation_y, k)
        end
    end

    scatter(bifurcation_x, bifurcation_y, title="Bifurcation Diagram", xlabel="p_jump", ylabel="Mean Oxidation State")
    savefig(joinpath(output_dir, "bifurcation_diagram.png"))
end

# Simulate System
function simulate(r::Int, initial_proteoform::String, steps::Int, p_jump::Float64)
    state, proteoforms = initialize_state(r, initial_proteoform)
    history = []
    entropies = []
    mean_oxidation_states = []

    for _ in 1:steps
        push!(history, deepcopy(state))
        entropy, mean_k = calc_metrics(state, proteoforms, r)
        push!(entropies, entropy)
        push!(mean_oxidation_states, mean_k)
        state = evolve_state(state, proteoforms, r, p_jump)
    end

    return history, proteoforms, entropies, mean_oxidation_states
end

# Main Execution
r = 3
initial_proteoform = "000"
steps = 1000
p_jump = 0.3
epsilon = 1e-5

# Run Simulation
history, proteoforms, entropies, mean_oxidation_states = simulate(r, initial_proteoform, steps, p_jump)
lyapunov = compute_lyapunov_exponent(r, initial_proteoform, steps, p_jump, epsilon)

println("Computed Lyapunov Exponent: ", lyapunov)

# Save Results
save_history(history, proteoforms)
save_metrics(steps, entropies, mean_oxidation_states)

# Generate Plots
plot_poincare(history, proteoforms, 10)
plot_bifurcation(r, initial_proteoform, steps, 0.1, 0.9, 50)
