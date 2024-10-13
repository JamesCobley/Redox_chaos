using XLSX
using Random

# Manually set general POX and PRED values (applies to all proteoforms)
POX = 0.7  # Probability of oxidation
PRED = 0.1  # Probability of reduction

# Load proteoform matrix from Excel file
function load_proteoforms(file_path::String)
    data = XLSX.readxlsx(file_path)
    sheet = data["Sheet1"]
    
    proteoforms = Dict()
    for row in 2:size(sheet, 1)  # Assuming row 1 has headers
        pf = sheet[row, 1]  # Proteoform ID
        allowed = split(sheet[row, 5], ", ")  # Allowed transitions
        barred = split(sheet[row, 6], ", ")  # Barred transitions
        proteoforms[pf] = (allowed, barred)
    end

    return proteoforms
end

# Function to run the Monte Carlo simulation
function run_simulation(proteoforms::Dict, molecules::Int, steps::Int)
    # Initialize population: All molecules start in PF001
    population = Dict{String, Int}()
    for _ in 1:molecules
        population["PF001"] = get(population, "PF001", 0) + 1
    end

    # Iterate over steps
    for step in 1:steps
        new_population = deepcopy(population)

        # Loop over the population and probabilistically apply transitions
        for (pf, count) in population
            allowed, _ = proteoforms[pf]  # Get allowed transitions for current proteoform

            for _ in 1:count
                rand_val = rand()

                if rand_val < POX  # Oxidation event
                    # Transition to a higher oxidation state
                    if !isempty(allowed)
                        new_pf = allowed[rand(1:length(allowed))]
                        new_population[pf] -= 1
                        new_population[new_pf] = get(new_population, new_pf, 0) + 1
                    end
                elseif rand_val < POX + PRED  # Reduction event
                    # Transition to a lower oxidation state
                    if !isempty(allowed)
                        new_pf = allowed[rand(1:length(allowed))]
                        new_population[pf] -= 1
                        new_population[new_pf] = get(new_population, new_pf, 0) + 1
                    end
                end
            end
        end

        population = new_population
    end

    return population
end

# Main function
function main()
    file_path = "path_to_your_excel_file.xlsx"  # Adjust this to your actual file path
    proteoforms = load_proteoforms(file_path)

    molecules = 100  # Number of molecules
    steps = 20  # Number of steps

    final_population = run_simulation(proteoforms, molecules, steps)
    
    println("Final Population Distribution:")
    for (pf, count) in final_population
        println("$pf: $count")
    end
end

main()
