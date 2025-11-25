# TorelliTrees

This [Julia](https://julialang.org/) package is used for the computations in _[paper title]_ by _[paper authors]_.
It produces [admcycles](https://pypi.org/project/admcycles/) code that, when executed, computes

$$\mathrm{Tor}^*[\mathcal{A}_{g_1}\times\cdots\times \mathcal{A}_{g_k}]\in R(\mathcal{M}_{g}^\mathrm{ct})$$

where $g=g_1+\dots+g_k$ and $\mathrm{Tor}:\mathcal{M}_g^\mathrm{ct}\rightarrow \mathcal{A}_G$ is the Torelli map.
It also illustrated the decorated stable trees indexing the relevant strata.


## Installation

To install the latest version of this package, first install [Julia](https://julialang.org/).
To install the `TorelliTrees` package, launch Julia and execute

    using Pkg
    Pkg.add(url="https://github.com/dh604/TorelliTrees.jl")

Once the package is installed, you can load it with

    using TorelliTrees

## Functionality

This package exports three functions: `stratum_trees`, `compute_contributions`, and `adm_code`.

### Step 1: Generating stratum trees using `stratum_trees`

This function returns the collection of stratum trees that are relevant for the computation of

$$\mathrm{Tor}^*[\mathcal{A}_{g_1}\times\cdots\times\mathcal{A}_{g_k}].$$


It also computes the smoothing relations between all of those trees.
For example, to generate the trees for $k=2, g_1=1, g_2=5$, use the following code in Julia:

    T15 = stratum_trees([1, 5])

#### Drawing pictures:
To export pdf files of the generated trees, use the optional argument `draw` together with the argument `folder_name`. This will produce LaTeX code to include the generated pictures into a LaTeX document, assuming the generated pictures will be placed in the folder `folder_name`.
For example, the code

    T15 = stratum_trees([1, 5]; draw=true, folder_name="A1A5/")

produces the relevant trees for $\mathrm{Tor}^*[\mathcal{A}_1\times\mathcal{A}_5]$ and produces pdf files `tree_1.pdf`, `tree_2.pdf`, etc. in the current directory.
It also prints the LaTeX code, assuming that those pdf files will be placed in the folder `A1A5/` relative to the LaTeX file.

### Step 2: Computing contributions using `compute_contributions`

To compute the polynomial $\mathrm{Cont}_T$ in the edge variables and Chern classes of the normal bundle for each $T$, simply run

    compute_contributions(T15)

### Step 3: Exporting `admcycles` code using `adm_code`

To obtain the admcycles code, execute

    adm_code(T15, "A1A5", "A1A5")

The first `"A1A5"` argument is a prefix used to name the variables in the sage file.
The second `"A1A5"` argument means that the output will be a file called `A1A5.sage`.

Executing this sage file will produce a variable called `Torelli_pullback`, which is precisely

$$\mathrm{Tor}^*[\mathcal{A}_{g_1}\times\cdots\times \mathcal{A}_{g_k}].$$

## Exampes

The folder `output` contains the resulting sage files for a number of small cases.