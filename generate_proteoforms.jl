# Function to generate binary proteoforms based on the number of cysteines (r)
function generate_proteoforms(r::Int)
    num_states = 2^r  # Total number of proteoforms
    proteoforms = [bitstring(i) for i in 0:(num_states - 1)]  # Generate binary strings for all states
    proteoforms = [lpad(p, r, '0')[end-r+1:end] for p in proteoforms]  # Ensure all strings are length r
    return proteoforms
end

# Function to display the generated proteoforms
function display_proteoforms(r::Int)
    proteoforms = generate_proteoforms(r)
    println("Generated proteoforms for r = $r:")
    for pf in proteoforms
        println(pf)
    end
end

# Prompt the user for the number of cysteines (r)
println("Enter the number of cysteines (r):")
r = parse(Int, readline())

# Display the proteoforms for the given r
display_proteoforms(r)
