@doc raw"""
    save_stratum_trees(filename::String, trees::Vector{StratumTree})

Save a vector of StratumTree objects to a file using Julia's serialization.

# Arguments
- `filename::String`: Path to the output file (will add .jls extension if not present)
- `trees::Vector{StratumTree}`: Vector of StratumTree objects to save

# Example
```julia
T14 = stratum_trees([1, 4])
save_stratum_trees("T14.jls", T14)
```

# Notes
The file format (.jls) is Julia-specific and preserves the exact structure of the objects,
including circular references. Files can only be loaded with compatible Julia and package versions.
"""
function save_stratum_trees(filename::String, STB::ST_backup)
    # Add .jls extension if not present
    if !endswith(filename, ".jls")
        filename = filename * ".jls"
    end

    open(filename, "w") do io
        serialize(io, STB)
        close(io)
    end
    
    trees = STB.STs
    cont_max = STB.last_index_with_cont_t
    println("Saved $(length(trees)) StratumTree(s) with $cont_max computed contributions to $filename")
    return filename
end

@doc raw"""
    load_stratum_trees(filename::String)

Load a vector of StratumTree objects from a file.

# Arguments
- `filename::String`: Path to the input file (will add .jls extension if not present)

# Returns
- `Vector{StratumTree}`: Vector of loaded StratumTree objects

# Example
```julia
T14 = load_stratum_trees("T14.jls")
```

# Notes
This function loads files created with `save_stratum_trees`. The file format is Julia-specific
and requires compatible Julia and package versions.
"""
function load_stratum_trees(filename::String)::ST_backup
    # Add .jls extension if not present
    if !endswith(filename, ".jls")
        filename = filename * ".jls"
    end

    STB = open(filename, "r") do io
        res = deserialize(io)
        close(io)
        return res
    end

    trees = STB.STs
    n_cont = STB.last_index_with_cont_t
    println("Loaded $(length(trees)) StratumTree(s) with $n_cont computed contributions from $filename")
    return STB
end
