# Redox Proteoform Chaos Simulation with Evolving P-Matrices
# Author: [James Cobley]
# Description: Simulation with 10 evolving P-matrices for distinct sub-populations,
# explicitly setting barred transitions to zero.

using Plots, Random, CSV, DataFrames, Statistics

# Set working directory
output_dir = pwd()

# Initialize system state with molecules
function initialize_state(r::Int, initial_proteoform::String, num_molecules::Int)
    proteoforms = [bitstring(i) |> s -> lpad(s, r, '0')[end-r+1:end] for i in 0:(2^r-1)]
    state = Dict(pf => 0.0 for pf in proteoforms)
    state[initial_proteoform] = num_molecules  # Start with specified number of molecules
    return state, proteoforms
end

# Create a P-matrix with barred transitions set to 0
function create_random_P_matrix(proteoforms::Vector{String})
    num_states = length(proteoforms)
    P = zeros(Float64, num_states, num_states)

    # Assign random probabilities to allowed transitions only
    for i in 1:num_states
        allowed_transitions = find_allowed_transitions(proteoforms[i], proteoforms, length(proteoforms[1]))
        
        # Assign random probabilities only to allowed transitions
        for pf in allowed_transitions
            j = findfirst(x -> x == pf, proteoforms)
            P[i, j] = rand()
        end

        # Normalize rows to ensure transition sum = 1
        row_sum = sum(P[i, :])
        if row_sum > 0
            P[i, :] /= row_sum
        end
    end

    return P
end

# Evolve multiple P-matrices for sub-populations
function evolve_multiple_P_matrices(
    state::Dict{String, Float64}, proteoforms::Vector{String}, P_matrices::Vector{Matrix{Float64}}
)
    new_state = Dict(pf => 0.0 for pf in proteoforms)
    sub_pool_size = sum(values(state)) / length(P_matrices)

    for (idx, P) in enumerate(P_matrices)
        for i in 1:length(proteoforms)
            for j in 1:length(proteoforms)
                new_state[proteoforms[j]] += state[proteoforms[i]] * P[i, j] * (state[proteoforms[i]] / sub_pool_size)
            end
        end
    end
    return new_state
end

# Calculate system metrics
function calc_metrics(state::Dict{String, Float64}, proteoforms::Vector{String}, r::Int)
    total_molecules = sum(values(state))
    probabilities = [state[pf] / total_molecules for pf in proteoforms]
    entropy = -sum(p * log(p + 1e-10) for p in probabilities)
    mean_k = sum(count('1', pf) * state[pf] for pf in proteoforms) / total_molecules
    return entropy, mean_k
end

# Save simulation data
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

# Generate Poincaré Plot
function plot_poincare(history::Vector{Dict{String, Float64}}, proteoforms::Vector{String}, period::Int)
    points_x = []
    points_y = []
    for t in 1:period:length(history)-period
        x = sum(count('1', pf) * history[t][pf] for pf in proteoforms)
        y = sum(count('1', pf) * history[t+period][pf] for pf in proteoforms)
        push!(points_x, x)
        push!(points_y, y)
    end
    scatter(points_x, points_y, title="Poincaré Diagram", xlabel="k(t)", ylabel="k(t+T)", xlim=(0, 8), ylim=(0, 8))
    savefig(joinpath(output_dir, "poincare_diagram.png"))
end

# Generate Bifurcation Plot
function plot_bifurcation(r::Int, initial_proteoform::String, steps::Int, p_min::Float64, p_max::Float64, p_steps::Int)
    p_values = range(p_min, p_max, length=p_steps)
    bifurcation_x = []
    bifurcation_y = []

    for p_jump in p_values
        _, proteoforms, _, mean_k_values = simulate_with_evolving_P_matrices(r, initial_proteoform, steps, 10_000)
        for k in mean_k_values[end-100:end]
            push!(bifurcation_x, p_jump)
            push!(bifurcation_y, k)
        end
    end

    scatter(bifurcation_x, bifurcation_y, title="Bifurcation Diagram", xlabel="p_jump", ylabel="Mean Oxidation State", xlim=(0.1, 0.9), ylim=(0, 8))
    savefig(joinpath(output_dir, "bifurcation_diagram.png"))
end

# Simulate system with evolving P-matrices
function simulate_with_evolving_P_matrices(r::Int, initial_proteoform::String, steps::Int, num_molecules::Int)
    state, proteoforms = initialize_state(r, initial_proteoform, num_molecules)
    num_pools = 10

    # Initialize 10 random P-matrices
    P_matrices = [create_random_P_matrix(proteoforms) for _ in 1:num_pools]

    history = Vector{Dict{String, Float64}}()
    entropies = Float64[]
    mean_oxidation_states = Float64[]

    for t in 1:steps
        # Update P-matrices every 100 steps
        if t % 100 == 0
            P_matrices = [create_random_P_matrix(proteoforms) for _ in 1:num_pools]
        end

        # Store system state
        push!(history, deepcopy(state))
        
        # Calculate metrics
        entropy, mean_k = calc_metrics(state, proteoforms, r)
        push!(entropies, entropy)
        push!(mean_oxidation_states, mean_k)

        # Evolve the system
        state = evolve_multiple_P_matrices(state, proteoforms, P_matrices)
    end

    return history, proteoforms, entropies, mean_oxidation_states
end

# Main Execution
r = 3
initial_proteoform = "000"
steps = 10_000
num_molecules = 10_000

history, proteoforms, entropies, mean_oxidation_states = simulate_with_evolving_P_matrices(
    r, initial_proteoform, steps, num_molecules
)

# Save Results
save_history(history, proteoforms)
save_metrics(steps, entropies, mean_oxidation_states)

# Generate Plots
plot_poincare(history, proteoforms, 10)
plot_bifurcation(r, initial_proteoform, steps, 0.1, 0.9, 50)
