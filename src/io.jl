@doc raw"""
    save_stratum_trees(filename::String, STB::ST_backup)

Save an ST_backup object to a file using Oscar's serialization.

# Arguments
- `filename::String`: Path to the output file (will add .oscar extension if not present)
- `STB::ST_backup`: ST_backup object to save

# Example
```julia
T14 = stratum_trees([1, 4])
save_stratum_trees("T14.oscar", STB)
```

# Notes
Uses Oscar's native serialization which properly handles polynomial types across
different Julia sessions and environments.
"""
function save_stratum_trees(filename::String, STB::ST_backup)
    # Add .oscar extension if not present
    if !endswith(filename, ".oscar")
        filename = filename * ".oscar"
    end

    Oscar.save(filename, STB)

    trees = STB.STs
    cont_max = STB.last_index_with_cont_t
    println("Saved $(length(trees)) StratumTree(s) with $cont_max computed contributions to $filename")
    return filename
end

@doc raw"""
    load_stratum_trees(filename::String)

Load an ST_backup object from a file.

# Arguments
- `filename::String`: Path to the input file (will add .oscar extension if not present)

# Returns
- `ST_backup`: Loaded ST_backup object

# Example
```julia
STB = load_stratum_trees("T14.oscar")
```

# Notes
This function loads files created with `save_stratum_trees`. Uses Oscar's native
serialization which properly handles polynomial types.
"""
function load_stratum_trees(filename::String)::ST_backup
    # Add .oscar extension if not present
    if !endswith(filename, ".oscar")
        filename = filename * ".oscar"
    end

    STB = Oscar.load(filename)

    trees = STB.STs
    n_cont = STB.last_index_with_cont_t
    println("Loaded $(length(trees)) StratumTree(s) with $n_cont computed contributions from $filename")
    return STB
end
