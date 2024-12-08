# Redox Proteoform Chaos Simulation with Evolving P-Matrices
# Author: [James Cobley]
# Description: Simulation with 10 evolving P-matrices for distinct sub-populations,
# explicitly setting barred transitions to zero.

using Plots, Random, CSV, DataFrames, Statistics

# Set working directory
output_dir = pwd()

# Initialize system state
function initialize_state(r::Int, initial_proteoform::String, num_molecules::Int)
    proteoforms = [bitstring(i) |> s -> lpad(s, r, '0')[end-r+1:end] for i in 0:(2^r-1)]
    state = Dict(pf => 0.0 for pf in proteoforms)
    state[initial_proteoform] = num_molecules
    return state, proteoforms
end

# Find allowed transitions
function find_allowed_transitions(proteoform::String, proteoforms::Vector{String}, r::Int)
    allowed = [proteoform]
    struct_vec = collect(proteoform)
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

# Create a P-matrix with safeguards
function create_random_P_matrix(proteoforms::Vector{String})
    num_states = length(proteoforms)
    P = zeros(Float64, num_states, num_states)

    for i in 1:num_states
        allowed_transitions = find_allowed_transitions(proteoforms[i], proteoforms, length(proteoforms[1]))
        for pf in allowed_transitions
            j = findfirst(x -> x == pf, proteoforms)
            P[i, j] = rand()
        end

        # Normalize rows with safeguard
        row_sum = sum(P[i, :])
        if row_sum > 0
            P[i, :] /= row_sum
        else
            println("Warning: Zero-row detected in P-matrix at row $i. Fixing to uniform.")
            P[i, :] .= 1.0 / num_states
        end
    end
    return P
end

# Evolve system using multiple P-matrices
function evolve_multiple_P_matrices(state::Dict{String, Float64}, proteoforms::Vector{String}, P_matrices::Vector{Matrix{Float64}})
    new_state = Dict(pf => 0.0 for pf in proteoforms)
    sub_pool_size = max(sum(values(state)), 1e-10) / length(P_matrices)

    for (idx, P) in enumerate(P_matrices)
        for i in 1:length(proteoforms)
            for j in 1:length(proteoforms)
                new_state[proteoforms[j]] += state[proteoforms[i]] * P[i, j] * (state[proteoforms[i]] / sub_pool_size)
            end
        end
    end

    # Normalize new state
    total_molecules = max(sum(values(new_state)), 1e-10)
    for pf in proteoforms
        new_state[pf] /= total_molecules
    end
    return new_state
end

# Simulate system with evolving P-matrices
function simulate_with_evolving_P_matrices(r::Int, initial_proteoform::String, steps::Int, num_molecules::Int)
    state, proteoforms = initialize_state(r, initial_proteoform, num_molecules)
    num_pools = 10
    P_matrices = [create_random_P_matrix(proteoforms) for _ in 1:num_pools]

    history = Vector{Dict{String, Float64}}()
    entropies = Float64[]
    mean_oxidation_states = Float64[]

    for t in 1:steps
        if t % 100 == 0
            P_matrices = [create_random_P_matrix(proteoforms) for _ in 1:num_pools]
        end

        push!(history, deepcopy(state))
        entropy, mean_k = calc_metrics(state, proteoforms, r)
        push!(entropies, entropy)
        push!(mean_oxidation_states, mean_k)

        state = evolve_multiple_P_matrices(state, proteoforms, P_matrices)
    end

    return history, proteoforms, entropies, mean_oxidation_states
end

# Calculate system metrics
function calc_metrics(state::Dict{String, Float64}, proteoforms::Vector{String}, r::Int)
    total_molecules = sum(values(state))
    probabilities = [state[pf] / total_molecules for pf in proteoforms]
    entropy = -sum(p * log(p + 1e-10) for p in probabilities)
    mean_k = sum(count('1', pf) * state[pf] for pf in proteoforms) / total_molecules
    return entropy, mean_k
end

# Compute Lyapunov Exponent
function compute_lyapunov_exponent(r::Int, initial_proteoform::String, steps::Int, num_molecules::Int, epsilon::Float64)
    original_state, proteoforms = initialize_state(r, initial_proteoform, num_molecules)
    perturbed_state = deepcopy(original_state)

    perturbed_state[proteoforms[2]] = min(perturbed_state[proteoforms[2]] + epsilon, num_molecules)
    perturbed_state[proteoforms[1]] = max(perturbed_state[proteoforms[1]] - epsilon, 0.0)

    P_matrices = [create_random_P_matrix(proteoforms) for _ in 1:10]
    distances = []

    for step in 1:steps
        dist = sqrt(sum((original_state[pf] - perturbed_state[pf])^2 for pf in proteoforms))
        push!(distances, max(dist, 1e-10))

        if step % 100 == 0
            P_matrices = [create_random_P_matrix(proteoforms) for _ in 1:10]
        end

        original_state = evolve_multiple_P_matrices(original_state, proteoforms, P_matrices)
        perturbed_state = evolve_multiple_P_matrices(perturbed_state, proteoforms, P_matrices)

        perturbed_total = max(sum(values(perturbed_state)), 1e-10)
        for pf in proteoforms
            perturbed_state[pf] /= perturbed_total
        end
    end

    return mean(log.(distances))
end

# Visualization Functions
function plot_final_state_distribution(state::Dict{String, Float64}, proteoforms::Vector{String})
    oxidation_levels = [count('1', pf) for pf in proteoforms]
    counts = [state[pf] for pf in proteoforms]
    bar(oxidation_levels, counts, xlabel="Oxidation Level (k)", ylabel="Molecule Count", title="Final State Distribution", legend=false)
    savefig("final_state_distribution.png")
end

function plot_poincare(history::Vector{Dict{String, Float64}}, proteoforms::Vector{String}, period::Int)
    points_x, points_y = [], []
    for t in 1:period:(length(history)-period)
        x = sum(count('1', pf) * history[t][pf] for pf in proteoforms)
        y = sum(count('1', pf) * history[t+period][pf] for pf in proteoforms)
        push!(points_x, x)
        push!(points_y, y)
    end
    scatter(points_x, points_y, title="Poincaré Diagram", xlabel="k(t)", ylabel="k(t+T)", legend=false)
    savefig("poincare_diagram.png")
end

# Main Execution
r = 3
initial_proteoform = "000"
steps = 10_000
num_molecules = 10_000
epsilon = 1e-5

# Run the Simulation
history, proteoforms, entropies, mean_oxidation_states = simulate_with_evolving_P_matrices(
    r, initial_proteoform, steps, num_molecules
)

# Compute and Print Lyapunov Exponent
lyapunov = compute_lyapunov_exponent(r, initial_proteoform, steps, num_molecules, epsilon)
println("Computed Lyapunov Exponent: ", lyapunov)

# Generate Plots
plot_final_state_distribution(history[end], proteoforms)
plot_poincare(history, proteoforms, 10)
