using XLSX  # For saving the output to Excel
using DataFrames  # For creating a DataFrame

# Function to generate binary proteoforms based on the number of cysteines (r)
function generate_proteoforms(r::Int)
    num_states = 2^r  # Total number of proteoforms
    proteoforms = [bitstring(i) for i in 0:(num_states - 1)]  # Generate binary strings for all states
    proteoforms = [lpad(p, r, '0')[end-r+1:end] for p in proteoforms]  # Ensure all strings are length r
    return proteoforms
end

# Function to find allowed transitions for each proteoform (including the "self" transition)
function find_allowed_transitions(proteoform::String, proteoforms::Vector{String}, r::Int)
    allowed = [proteoform]  # Include the proteoform itself as an allowed transition (self)
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
# Barred transitions are all proteoforms not in the allowed set
function find_barred_transitions(proteoform::String, allowed_transitions::Vector{String}, proteoforms::Vector{String})
    return [pf for pf in proteoforms if !(pf in allowed_transitions)]
end

# Function to save the data into an Excel file
function save_to_excel(r::Int, data::DataFrame, file_name::String)
    # Open or create the Excel file
    XLSX.openxlsx(file_name, mode="w") do file
        # Write the DataFrame contents to the sheet
        sheet_name = "Proteoforms"
        XLSX.writetable(file, collect(eachcol(data)), header=names(data), sheetname=sheet_name)
    end
    println("Excel file saved successfully as $file_name.")
end

# Function to display and save the proteoforms, allowed/barred transitions, and conservation of degrees
function display_and_save_transitions(r::Int, file_name::String)
    proteoforms = generate_proteoforms(r)
    println("Generated proteoforms for r = $r:")
    
    # Initialize a DataFrame to store the data
    data = DataFrame(PF = String[], k_value = Int[], Percent_OX = Float64[], Structure = String[], 
                     Allowed = String[], Barred = String[], K_minus_0 = Int[], K_plus = Int[], Conservation_of_degrees = Int[])

    for (i, proteoform) in enumerate(proteoforms)
        k_value = count(c -> c == '1', proteoform)  # Count the number of '1's (oxidized cysteines)
        percent_ox = 100 * k_value / r  # Calculate the percentage of oxidation

        # Find allowed transitions
        allowed_transitions = find_allowed_transitions(proteoform, proteoforms, r)
        barred_transitions = find_barred_transitions(proteoform, allowed_transitions, proteoforms)

        # Calculate Conservation of Degrees (r + 1 includes self)
        K_minus_0 = sum(k < k_value for k in [count(c -> c == '1', pf) for pf in allowed_transitions])
        K_plus = sum(k > k_value for k in [count(c -> c == '1', pf) for pf in allowed_transitions])
        conservation_of_degrees = K_minus_0 + K_plus + 1  # Adding 1 for the "self" transition

        # Append the results to the DataFrame
        push!(data, (PF = "PF$(lpad(i, 3, '0'))", k_value = k_value, Percent_OX = percent_ox, 
                     Structure = proteoform, Allowed = join(allowed_transitions, ", "), 
                     Barred = join(barred_transitions, ", "), K_minus_0 = K_minus_0, 
                     K_plus = K_plus, Conservation_of_degrees = conservation_of_degrees))
    end

    # Display the generated DataFrame
    println(data)

    # Save the DataFrame to an Excel file
    save_to_excel(r, data, file_name)
end

# Prompt the user for the number of cysteines (r)
println("Enter the number of cysteines (r):")
r = parse(Int, readline())

# Prompt the user for the file name to save the Excel file
println("Enter the file name (with .xlsx extension) to save the output:")
file_name = readline()

# Display and save the proteoforms and transitions for the given r
display_and_save_transitions(r, file_name)
