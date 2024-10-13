using DataFrames
using XLSX  # For saving the output as an Excel file

# Function to generate the proteoform dictionary based on the number of cysteines (r)
function generate_proteoform_dict(r::Int)
    num_states = 2^r  # Total possible proteoforms (2^r)
    proteoforms = [bitstring(i) for i in 0:(num_states - 1)]  # Generate binary strings for all states
    proteoforms = [lpad(p, r, '0') for p in proteoforms]  # Ensure all strings are length r (capped at r cysteines)

    # Create dictionary: map PF IDs to binary structures
    proteoform_dict = Dict{String, String}()  # Explicitly typed as Dict{String, String}
    for i in 1:num_states
        PF = "PF$(lpad(i, 3, '0'))"  # PF001, PF002, ...
        structure = join([string(c) for c in split(proteoforms[i], "")], ",")  # Binary structure like "0,0", "1,0", etc.
        proteoform_dict[PF] = structure
    end
    return proteoform_dict
end

# Function to generate allowed and barred transitions
function generate_transitions(proteoform_dict::Dict{String, String}, r::Int)
    PF = collect(keys(proteoform_dict))  # List of proteoform IDs (e.g., PF001, PF002)
    num_states = length(PF)
    k_value = [count(c -> c == '1', proteoform_dict[pf]) for pf in PF]  # Count of oxidized cysteines (1s)
    percent_ox = [100 * k / r for k in k_value]  # Percentage of oxidation

    # Function to find allowed transitions (stepwise +1 or -1 in k_value, considering both positions)
    function find_allowed_transitions(current_structure::String)
        current_k = count(c -> c == '1', split(current_structure, ","))  # Count current oxidized cysteines
        allowed = []
        struct_vec = split(current_structure, ",")  # Turn structure into an array of substrings

        # Check for stepwise changes at each cysteine position
        for i in 1:r
            # Toggle the i-th cysteine between oxidized ("1") and reduced ("0")
            new_structure = vcat(struct_vec[1:i-1], (struct_vec[i] == "0" ? "1" : "0"), struct_vec[i+1:end])
            new_k = count(c -> c == '1', new_structure)  # Recompute the oxidation count

            # Only allow transitions with ±1 change in k_value and no side swaps
            if abs(new_k - current_k) == 1
                new_struct_str = join(new_structure, ",")  # Convert back to "0,1", etc.
                allowed_pf = findfirst(x -> proteoform_dict[x] == new_struct_str, PF)
                if allowed_pf !== nothing  # Only push valid transitions
                    push!(allowed, PF[allowed_pf])
                end
            end
        end
        return allowed
    end

    # Generate allowed and barred transitions
    allowed = [join(find_allowed_transitions(proteoform_dict[pf]), ", ") for pf in PF]

    # Function to find barred transitions
    function find_barred_transitions(current_pf::String, allowed_pfs::Vector{String})
        return [pf for pf in PF if !(pf in allowed_pfs) && pf != current_pf]
    end

    barred = [join(find_barred_transitions(PF[i], [String(x) for x in split(allowed[i], ", ")]), ", ") for i in 1:num_states]

    # Calculate K - 0 and K + for transitions
    K_minus_0 = [sum(k < k_value[i] for k in [count(c -> c == '1', proteoform_dict[pf]) for pf in split(allowed[i], ", ")]) for i in 1:num_states]
    K_plus = [sum(k > k_value[i] for k in [count(c -> c == '1', proteoform_dict[pf]) for pf in split(allowed[i], ", ")]) for i in 1:num_states]

    # Ensure conservation of degrees is equal to r
    conservation_of_degrees = [K_minus_0[i] + K_plus[i] for i in 1:num_states]

    # Create the DataFrame
    df = DataFrame(
        "PF" => PF,
        "k_value" => k_value,
        "Percent (OX)" => percent_ox,
        "Structure" => collect(values(proteoform_dict)),  # Proteoform structure from the dictionary
        "Allowed" => allowed,
        "Barred" => barred,
        "K - 0" => K_minus_0,
        "K +" => K_plus,
        "Conservation of degrees" => conservation_of_degrees
    )
    
    return df
end

# Function to save DataFrame to Excel
function save_to_excel(df::DataFrame, file_path::String)
    try
        XLSX.openxlsx(file_path, mode="w") do workbook
            sheet = XLSX.addsheet(workbook, "Proteoforms")
            XLSX.writetable(sheet, collect(eachcol(df)), header=names(df))
        end
        println("Excel file saved as $file_path")
    catch e
        println("Error saving the file: $e")
    end
end

# Function to run the script in a terminal without widgets
function run_in_terminal()
    # Prompt the user for the number of cysteines (r)
    println("Enter the number of cysteines (r):")
    r = parse(Int, readline())

    # Generate the proteoform dictionary
    proteoform_dict = generate_proteoform_dict(r)

    # Generate transitions and build the DataFrame
    df = generate_transitions(proteoform_dict, r)

    # Prompt the user for the file path to save the output
    println("Enter the file path where you want to save the output:")
    file_path = readline()

    # Ensure that the path includes the filename
    if !endswith(file_path, ".xlsx")
        file_path *= "/proteoform_transitions_r_$(r).xlsx"
    end

    # Display the DataFrame
    println("Generated DataFrame:")
    display(df)

    # Save the DataFrame to an Excel file
    save_to_excel(df, file_path)
end

# Call the function to run in terminal
run_in_terminal()
