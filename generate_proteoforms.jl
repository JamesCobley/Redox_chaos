using DataFrames
using Interact  # For creating the interactive widget
using XLSX      # For saving the output as an Excel file

# Function to generate the proteoform transitions based on the number of cysteines (r)
function generate_proteoforms(r::Int)
    # Generate all possible proteoforms (2^r combinations)
    num_states = 2^r
    proteoforms = [bitstring(i) for i in 0:(num_states - 1)]
    proteoforms = [lpad(p, r, '0') for p in proteoforms]  # Ensure all strings are length r

    # Define DataFrame columns
    PF = ["PF$(lpad(i + 1, 3, '0'))" for i in 0:(num_states - 1)]
    k_value = [count(c -> c == '1', p) for p in proteoforms]  # Count of oxidized cysteines (1s)
    percent_ox = [100 * k / r for k in k_value]  # Percentage of oxidation

    # Binary structure of proteoforms (0=reduced, 1=oxidized)
    structure = [replace(p, '0' => '0', '1' => '1') for p in proteoforms]

    # Function to find allowed transitions (+/- 1 oxidation/reduction)
    function find_allowed_transitions(struct, r)
        current_k = count(c -> c == '1', struct)
        allowed = []
        # Check transitions where one site changes (either oxidation or reduction)
        for i in 1:r
            new_struct = struct[1:i-1] * (struct[i] == '0' ? "1" : "0") * struct[i+1:end]
            new_k = count(c -> c == '1', new_struct)
            if abs(new_k - current_k) == 1  # Only allow stepwise transitions
                push!(allowed, new_struct)
            end
        end
        return allowed
    end

    # Generate allowed and barred transitions
    allowed = [join(["PF$(lpad(parse(Int, b, base=2) + 1, 3, '0'))" for b in find_allowed_transitions(p, r)], ", ") for p in proteoforms]

    # Function to find barred transitions
    function find_barred_transitions(current_pf, allowed_pfs)
        return [pf for pf in PF if !(pf in allowed_pfs) && pf != current_pf]
    end

    barred = [join(find_barred_transitions(PF[i], split(allowed[i], ", ")), ", ") for i in 1:num_states]

    # Define columns for transitions
    K_minus_0 = [sum(k < current_k for k in [count(c -> c == '1', p) for p in split(allowed[i], ", ")]) for i in 1:num_states]
    K_plus = [sum(k > current_k for k in [count(c -> c == '1', p) for p in split(allowed[i], ", ")]) for i in 1:num_states]
    conservation_of_degrees = [K_minus_0[i] + K_plus[i] for i in 1:num_states]

    # Create the DataFrame
    df = DataFrame(
        PF = PF,
        k_value = k_value,
        "Percent (OX)" = percent_ox,
        Structure = structure,
        Allowed = allowed,
        Barred = barred,
        "K - 0" = K_minus_0,
        "K +" = K_plus,
        "Conservation of degrees" = conservation_of_degrees
    )
    
    return df
end

# Interactive widget using Interact.jl
@manipulate for r in 1:10, file_path = textbox("Enter the file path:")
    # Generate the proteoform transitions based on the selected r
    df = generate_proteoforms(r)

    # Display the DataFrame
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

