using DataFrames
using XLSX  # For saving the output as an Excel file

# Function to generate the proteoform dictionary based on the number of cysteines (r)
function generate_proteoform_dict(r::Int)
    num_states = 2^r  # Total possible proteoforms (2^r)
    proteoforms = [bitstring(i) for i in 0:(num_states - 1)]  # Generate binary strings for all states
    proteoforms = [lpad(p, r, '0')[end-r+1:end] for p in proteoforms]  # Ensure all strings are length r

    # Create dictionary: map PF IDs to binary structures
    proteoform_dict = Dict{String, String}()
    for i in 1:num_states
        PF = "PF$(lpad(i, 3, '0'))"  # Proteoform ID (e.g., PF001, PF002)
        structure = join([string(c) for c in split(proteoforms[i], "")], ",")  # Binary structure, separated by commas
        proteoform_dict[PF] = structure
    end
    return proteoform_dict
end

# Function to find allowed transitions and directly structure them with `;`
function find_allowed_transitions(proteoform_dict::Dict{String, String}, current_pf::String, r::Int)
    current_structure = proteoform_dict[current_pf]
    struct_vec = split(current_structure, ",")  # Split into array for each cysteine state
    allowed = []

    current_k = count(c -> c == "1", struct_vec)  # Number of oxidized cysteines

    # Loop over each cysteine and toggle between oxidized ("1") and reduced ("0")
    for i in 1:r
        new_structure = copy(struct_vec)  # Copy the structure
        new_structure[i] = (new_structure[i] == "0" ? "1" : "0")  # Toggle the i-th cysteine
        new_k = count(c -> c == "1", new_structure)  # Recalculate the oxidation count

        # Only allow stepwise transitions (±1 in k_value)
        if abs(new_k - current_k) == 1
            new_struct_str = join(new_structure, ",")

            # Lookup the proteoform with the new structure, ensuring keys exist
            allowed_pf = findfirst(x -> haskey(proteoform_dict, x) && proteoform_dict[x] == new_struct_str, collect(keys(proteoform_dict)))
            if allowed_pf !== nothing && allowed_pf != ""  # Skip invalid or non-existent keys
                push!(allowed, allowed_pf)
            end
        end
    end

    return join(allowed, ";")  # Return allowed transitions separated by `;`
end

# Function to find barred transitions and directly structure them with `;`
function find_barred_transitions(current_pf::String, allowed_pfs::String, PF::Vector{String})
    allowed_list = split(allowed_pfs, ";")  # Split allowed proteoforms by `;`
    return join([pf for pf in PF if !(pf in allowed_list) && pf != current_pf], ";")
end

# Function to generate allowed and barred transitions
function generate_transitions(proteoform_dict::Dict{String, String}, r::Int)
    PF = collect(keys(proteoform_dict))  # List of proteoform IDs (e.g., PF001, PF002)
    num_states = length(PF)
    k_value = [count(c -> c == '1', proteoform_dict[pf]) for pf in PF]  # Count of oxidized cysteines (1s)
    percent_ox = [100 * k / r for k in k_value]  # Percentage of oxidation

    # Generate allowed and barred transitions for each proteoform
    allowed = [find_allowed_transitions(proteoform_dict, PF[i], r) for i in 1:num_states]
    barred = [find_barred_transitions(PF[i], allowed[i], PF) for i in 1:num_states]

    # Calculate K - 0 and K + for transitions
    K_minus_0 = [sum(k < k_value[i] for k in [count(c -> c == '1', proteoform_dict[pf]) for pf in split(allowed[i], ";")]) for i in 1:num_states]
    K_plus = [sum(k > k_value[i] for k in [count(c -> c == '1', proteoform_dict[pf]) for pf in split(allowed[i], ";")]) for i in 1:num_states]

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
