# Input Scripts

In the following the expected schemes of the input scripts is described.
Therefore, each section of the input script gets its own chapter. All sections
should return a `list` as results, which are concatenated to one setup list at
the end of the file. All values must have the same name in order to get parsed
by POET.

## Grid initialization

| name           | type           | description                                                           |
|----------------|----------------|-----------------------------------------------------------------------|
| `n_cells`      | Numeric Vector | Number of cells in each direction                                     |
| `s_cells`      | Numeric Vector | Spatial resolution of grid in each direction                          |
| `type`         | String         | Type of initialization, can be set to *scratch*, *phreeqc* or *rds*   |

## Diffusion parameters

| name           | type                 | description                               |
|----------------|----------------------|-------------------------------------------|
| `init`         | Named Numeric Vector | Initial state for each diffused species   |
| `vecinj`       | Data Frame           | Defining all boundary conditions row wise |
| `vecinj_inner` | List of Triples      | Inner boundaries                          |
| `vecinj_index` | List of 4 elements   | Ghost nodes boundary conditions           |
| `alpha`        | Named Numeric Vector | Constant alpha for each species           |

### Remark on boundary conditions

Each boundary condition should be defined in `vecinj` as a data frame, where one
row holds one boundary condition.

To define inner (constant) boundary conditions, use a list of triples in
`vecinj_inner`, where each triples is defined by $(i,x,y)$. $i$ is defining the
boundary condition, referencing to the row in `vecinj`. $x$ and $y$ coordinates
then defining the position inside the grid. 

Ghost nodes are set by `vecinj_index` which is a list containing boundaries for
each celestial direction (**important**: named by `N, E, S, W`). Each direction
is a numeric vector, also representing a row index of the `vecinj` data frame
for each ghost node, starting at the left-most and upper cell respectively. By
setting the boundary condition to $0$, the ghost node is set as closed boundary.

#### Example

Suppose you have a `vecinj` data frame defining 2 boundary conditions and a grid
consisting of $10 \times 10$ grid cells. Grid cell $(1,1)$ should be set to the
first boundary condition and $(5,6)$ to the second. Also, all boundary
conditions for the ghost nodes should be closed. Except the southern boundary,
which should be set to the first boundary condition injection. The following
setup describes how to setup your initial script, where `n` and `m` are the
grids cell count for each direction ($n = m = 10$):

```R
vecinj_inner <- list (
  l1 = c(1, 1, 1),
  l2 = c(2, 5, 6)
)

vecinj_index <- list(
  "N" = rep(0, n),
  "E" = rep(0, m),
  "S" = rep(1, n),
  "W" = rep(0, m)
)
```

## Chemistry parameters

| name           | type         | description                                                                      |
|----------------|--------------|----------------------------------------------------------------------------------|
| `database`     | String       | Path to the Phreeqc database                                                     |
| `input_script` | String       | Path the the Phreeqc input script                                                |
| `dht_species`  | Named Vector | Indicates significant digits to use for each species for DHT rounding.           |
| `pht_species`  | Named Vector | Indicates significant digits to use for each species for Interpolation rounding. |

## Final setup

| name           | type           | description                                                |
|----------------|----------------|------------------------------------------------------------|
| `grid`         | List           | Grid parameter list                                        |
| `diffusion`    | List           | Diffusion parameter list                                   |
| `chemistry`    | List           | Chemistry parameter list                                   |
| `iterations`   | Numeric Value  | Count of iterations                                        |
| `timesteps`    | Numeric Vector | $\Delta t$ to use for specific iteration                   |
| `store_result` | Boolean        | Indicates if results should be stored                      |
| `out_save`     | Numeric Vector | *optional:* At which iteration the states should be stored |
