import "regent"

--TODO: #67: https://github.com/stfc/RegentParticleDSL/issues/67
--If a particle leaves the domain it is just treated as though its still in the domain. The physics implementation
--is expected to avoid this happening as it is undefined behaviour at this moment
if DSL_DIMENSIONALITY == 2 then
  neighbour_init = require("src/neighbour_search/cell_pair_tradequeues_nonperiod/2d_neighbour_init")
elseif DSL_DIMENSIONALITY == 3 then
  neighbour_init = require("src/neighbour_search/cell_pair_tradequeues_nonperiod/3d_neighbour_init")
end

return neighbour_init
