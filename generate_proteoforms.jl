using DataFrames
using XLSX  # For saving the output as an Excel file

# Function to generate the proteoform transitions based on the number of cysteines (r)
function generate_proteoforms(r::Int)
    # Generate all possible proteoforms (2^r combinations)
    num_states = 2^r
    proteoforms = [bitstring(i) for i in 0:(num_states - 1)]
    proteoforms = [lpad(p, r, '0') for p in proteoforms]  # Ensure all strings are length r

    # Define DataFrame columns
    PF = ["PF$(lpad(BigInt(parse(BigInt, p, base=2)) + 1, 3, '0'))" for p in proteoforms]  # Use BigInt here
    k_value = [count(c -> c == '1', p) for p in proteoforms]  # Count of oxidized cysteines (1s)
    percent_ox = [100 * k / r for k in k_value]  # Percentage of oxidation

    # Binary structure of proteoforms (0=reduced, 1=oxidized)
    structure = [replace(p, '0' => '0', '1' => '1') for p in proteoforms]

    # Function to find allowed transitions (based on stepwise oxidation/reduction)
    function find_allowed_transitions(structure::String, r::Int)
        current_k = count(c -> c == '1', structure)  # Count of current oxidized cysteines
        allowed = []
        # Check transitions where one site changes (either oxidation or reduction)
        for i in 1:r
            # Toggle the i-th cysteine between oxidized ("1") and reduced ("0")
            new_structure = structure[1:i-1] * (structure[i] == '0' ? "1" : "0") * structure[i+1:end]
            new_k = count(c -> c == '1', new_structure)  # Recompute the oxidation count
            if abs(new_k - current_k) == 1  # Only allow transitions where one oxidation/reduction occurs
                push!(allowed, new_structure)
            end
        end
        return allowed
    end

    # Generate allowed and barred transitions
    allowed = [join(["PF$(lpad(BigInt(parse(BigInt, b, base=2)) + 1, 3, '0'))" for b in find_allowed_transitions(p, r)], ", ") for p in proteoforms]

    # Function to find barred transitions
    function find_barred_transitions(current_pf::String, allowed_pfs::Vector{String})
        return [pf for pf in PF if !(pf in allowed_pfs) && pf != current_pf]
    end

    # Convert SubString to String in split results to avoid type mismatch
    barred = [join(find_barred_transitions(PF[i], [String(x) for x in split(allowed[i], ", ")]), ", ") for i in 1:num_states]

    # Define columns for transitions
    K_minus_0 = [sum(k < count(c -> c == '1', proteoforms[i]) for k in [count(c -> c == '1', p) for p in split(allowed[i], ", ")]) for i in 1:num_states]
    K_plus = [sum(k > count(c -> c == '1', proteoforms[i]) for k in [count(c -> c == '1', p) for p in split(allowed[i], ", ")]) for i in 1:num_states]
    conservation_of_degrees = [K_minus_0[i] + K_plus[i] for i in 1:num_states]

    # Create the DataFrame
    df = DataFrame(
        PF = PF,
        k_value = k_value,
        Percent_OX = percent_ox,  # Changed to valid symbol format
        Structure = structure,
        Allowed = allowed,
        Barred = barred,
        K_minus_0 = K_minus_0,
        K_plus = K_plus,
        Conservation_of_degrees = conservation_of_degrees
    )
    
    return df
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
    if !isempty(file_path)
        file_name = "$(file_path)/proteoform_transitions_r_$(r).xlsx"
        try
            XLSX.openxlsx(file_name, mode="w") do workbook
                sheet = XLSX.addsheet(workbook, "Proteoforms")
                XLSX.writetable(sheet, collect(eachcol(df)), header=names(df))
            end
            println("Excel file saved as $file_name")
        catch e
            println("Error saving the file: $e")
        end
    else
        println("Please provide a valid file path to save the file.")
    end
end

# Call the function to run in terminal
run_in_terminal()
