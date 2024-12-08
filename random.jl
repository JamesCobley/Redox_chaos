# Redox Proteoform Chaos Simulation with Evolving P-Matrices
# Author: [James Cobley]
# Description: Simulation with 10 evolving P-matrices for distinct sub-populations,
# explicitly setting barred transitions to zero.

# Corrected Full Code with Proper Ends and Lyapunov Exponent
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

# Create a P-matrix with barred transitions set to 0
function create_random_P_matrix(proteoforms::Vector{String})
    num_states = length(proteoforms)
    P = zeros(Float64, num_states, num_states)

    for i in 1:num_states
        allowed_transitions = find_allowed_transitions(proteoforms[i], proteoforms, length(proteoforms[1]))
        for pf in allowed_transitions
            j = findfirst(x -> x == pf, proteoforms)
            P[i, j] = rand()  # Assign random probabilities only to allowed transitions
        end

        # Normalize rows
        row_sum = sum(P[i, :])
        if row_sum > 0
            P[i, :] /= row_sum
        end
    end

    return P
end

# Evolve system using multiple P-matrices
function evolve_multiple_P_matrices(state::Dict{String, Float64}, proteoforms::Vector{String}, P_matrices::Vector{Matrix{Float64}})
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

# Compute the Lyapunov Exponent
function compute_lyapunov_exponent(r::Int, initial_proteoform::String, steps::Int, num_molecules::Int, epsilon::Float64)
    original_state, proteoforms = initialize_state(r, initial_proteoform, num_molecules)
    perturbed_state = deepcopy(original_state)
    perturbed_state[proteoforms[2]] += epsilon
    perturbed_state[proteoforms[1]] -= epsilon

    P_matrices = [create_random_P_matrix(proteoforms) for _ in 1:10]

    distances = []

    for step in 1:steps
        dist = sqrt(sum((original_state[pf] - perturbed_state[pf])^2 for pf in proteoforms))
        push!(distances, dist)

        # Update P-matrices every 100 steps
        if step % 100 == 0
            P_matrices = [create_random_P_matrix(proteoforms) for _ in 1:10]
        end

        original_state = evolve_multiple_P_matrices(original_state, proteoforms, P_matrices)
        perturbed_state = evolve_multiple_P_matrices(perturbed_state, proteoforms, P_matrices)

        # Normalize the perturbed state
        perturbed_total = sum(values(perturbed_state))
        for pf in proteoforms
            perturbed_state[pf] /= perturbed_total
        end
    end

    return mean(log.(distances .+ 1e-10))
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

