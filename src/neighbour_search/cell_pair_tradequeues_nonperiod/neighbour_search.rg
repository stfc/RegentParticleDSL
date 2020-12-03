import "regent"

if DSL_DIMENSIONALITY == 2 then
  require("src/neighbour_search/cell_pair_tradequeues_nonperiod/2d_neighbour_search")
else if DSL_DIMENSIONALITY == 3 then
  require("src/neighbour_search/cell_pair_tradequeues_nonperiod/3d_neighbour_search")
end
