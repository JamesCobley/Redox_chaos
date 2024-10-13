using DataFrames
using XLSX  # For saving the output as an Excel file

# Function to generate the proteoform transitions based on the number of cysteines (r)
function generate_proteoforms(r::Int)
    # Generate all possible proteoforms (2^r combinations)
    num_states = 2^r
    proteoforms = [bitstring(i) for i in 0:(num_states - 1)]
    proteoforms = [lpad(p, r, '0') for p in proteoforms]  # Ensure all strings are length r

    # Define DataFrame columns
    PF = ["PF$(lpad(i + 1, 3, '0'))" for i in 0:(num_states - 1)]  # Sequential numbering: PF001, PF002, etc.
    k_value = [count(c -> c == '1', p) for p in proteoforms]  # Count of oxidized cysteines (1s)
    percent_ox = [100 * k / r for k in k_value]  # Percentage of oxidation

    # Binary structure of proteoforms (0=reduced, 1=oxidized)
    structure = [join([string(c) for c in p], ",") for p in proteoforms]  # Turn into "0,0", "1,1", etc.

    # Function to find allowed transitions (based on stepwise oxidation/reduction)
    function find_allowed_transitions(structure::String, r::Int)
        current_k = count(c -> c == '1', split(structure, ","))  # Count of current oxidized cysteines
        allowed = []
        struct_vec = split(structure, ",")  # Turn the structure into an array of substrings

        # Check transitions where one site changes (either oxidation or reduction)
        for i in 1:r
            # Toggle the i-th cysteine between oxidized ("1") and reduced ("0")
            new_structure = vcat(struct_vec[1:i-1], (struct_vec[i] == "0" ? "1" : "0"), struct_vec[i+1:end])
            new_k = count(c -> c == '1', new_structure)  # Recompute the oxidation count
            if abs(new_k - current_k) == 1  # Only allow transitions where one oxidation/reduction occurs
                push!(allowed, join(new_structure, ","))  # Join new structure back into "0,1", etc.
            end
        end
        return allowed
    end

    # Generate allowed and barred transitions
    allowed = [join(["PF$(lpad(findfirst(x -> x == b, structure) + 1, 3, '0'))" for b in find_allowed_transitions(p, r)], ", ") for p in structure]

    # Function to find barred transitions
    function find_barred_transitions(current_pf::String, allowed_pfs::Vector{String})
        return [pf for pf in PF if !(pf in allowed_pfs) && pf != current_pf]
    end

    barred = [join(find_barred_transitions(PF[i], [String(x) for x in split(allowed[i], ", ")]), ", ") for i in 1:num_states]

    # Define columns for transitions
    K_minus_0 = [sum(k < count(c -> c == '1', proteoforms[i]) for k in [count(c -> c == '1', p) for p in split(allowed[i], ", ")]) for i in 1:num_states]
    K_plus = [sum(k > count(c -> c == '1', proteoforms[i]) for k in [count(c -> c == '1', p) for p in split(allowed[i], ", ")]) for i in 1:num_states]
    conservation_of_degrees = [K_minus_0[i] + K_plus[i] for i in 1:num_states]

    # Correct DataFrame creation, using the => operator for all columns
    df = DataFrame(
        "PF" => PF,
        "k_value" => k_value,
        "Percent (OX)" => percent_ox,  # Correct syntax for column with parentheses
        "Structure" => structure,
        "Allowed" => allowed,
        "Barred" => barred,
        "K - 0" => K_minus_0,          # Correct syntax for column with special characters
        "K +" => K_plus,               # Correct syntax for column with special characters
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

    # Prompt the user for the file path to save the output
    println("Enter the file path where you want to save the output:")
    file_path = readline()

    # Generate the proteoform transitions based on the selected r
    df = generate_proteoforms(r)

    # Display the DataFrame
    println("Generated DataFrame:")
    display(df)

    # Save the DataFrame to an Excel file
    save_to_excel(df, file_path)
end

# Call the function to run in terminal
run_in_terminal()
