# Function to generate binary proteoforms based on the number of cysteines (r)
function generate_proteoforms(r::Int)
    num_states = 2^r  # Total number of proteoforms
    proteoforms = [bitstring(i) for i in 0:(num_states - 1)]  # Generate binary strings for all states
    proteoforms = [lpad(p, r, '0')[end-r+1:end] for p in proteoforms]  # Ensure all strings are length r
    return proteoforms
end

# Function to find allowed transitions for each proteoform
function find_allowed_transitions(proteoform::String, proteoforms::Vector{String}, r::Int)
    allowed = []
    struct_vec = collect(proteoform)  # Split the proteoform into an array of '0's and '1's

    # Current k_value (number of oxidized cysteines)
    current_k = count(c -> c == '1', struct_vec)

    # Check transitions by toggling one cysteine at a time
    for i in 1:r
        new_structure = copy(struct_vec)
        new_structure[i] = (new_structure[i] == '0' ? '1' : '0')  # Toggle cysteine oxidation state
        new_proteoform = join(new_structure)
        new_k = count(c -> c == '1', new_structure)

        # Only allow transitions where k_value changes by ±1
        if abs(new_k - current_k) == 1
            push!(allowed, new_proteoform)
        end
    end

    return allowed
end

# Function to find barred transitions for each proteoform
function find_barred_transitions(proteoform::String, allowed_transitions::Vector{String}, proteoforms::Vector{String})
    return [pf for pf in proteoforms if !(pf in allowed_transitions) && pf != proteoform]
end

# Function to display the proteoforms, allowed/barred transitions, and conservation of degrees
function display_transitions(r::Int)
    proteoforms = generate_proteoforms(r)
    println("Generated proteoforms for r = $r:")
    
    for (i, proteoform) in enumerate(proteoforms)
        k_value = count(c -> c == '1', proteoform)  # Count the number of '1's (oxidized cysteines)
        percent_ox = 100 * k_value / r  # Calculate the percentage of oxidation

        # Find allowed transitions
        allowed_transitions = find_allowed_transitions(proteoform, proteoforms, r)
        barred_transitions = find_barred_transitions(proteoform, allowed_transitions, proteoforms)

        # Calculate Conservation of Degrees
        K_minus_0 = sum(k < k_value for k in [count(c -> c == '1', pf) for pf in allowed_transitions])
        K_plus = sum(k > k_value for k in [count(c -> c == '1', pf) for pf in allowed_transitions])
        conservation_of_degrees = K_minus_0 + K_plus

        # Print information for each proteoform
        println("PF$(lpad(i, 3, '0')): $proteoform, k_value: $k_value, Percent Oxidation: $percent_ox%")
        println("  Allowed transitions: $(join(allowed_transitions, ", "))")  # Fixed the join function
        println("  Barred transitions: $(join(barred_transitions, ", "))")   # Fixed the join function
        println("  K - 0: $K_minus_0, K +: $K_plus, Conservation of degrees: $conservation_of_degrees")
        println()
    end

    # Print the total number of proteoforms
    println("\nTotal number of proteoforms: $(length(proteoforms))")
end

# Prompt the user for the number of cysteines (r)
println("Enter the number of cysteines (r):")
r = parse(Int, readline())

# Display the proteoforms and transitions for the given r
display_transitions(r)

