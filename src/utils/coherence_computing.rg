-------------------------------------------------------------
--Copyright 2020 Science and Technologies Facilities Council
--Licensed under the MIT License
--Author Aidan Chalk, STFC Hartree Centre

import "regent"

local coherence_compute = {}

function coherence_compute.compute_coherences_pair_task(update_neighbours, parts1, parts2)
    local coherences = terralib.newlist()
    if update_neighbours then
      coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
      coherences:insert( regentlib.coherence( regentlib.exclusive, parts2 ) )
    else
      coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
      coherences:insert( regentlib.coherence( regentlib.exclusive, parts2 ) )
    end
    return coherences
end

function coherence_compute.compute_coherences_self_task(update_neighbours, parts1)
    local coherences = terralib.newlist()
    if update_neighbours then
      coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
    else
      coherences:insert( regentlib.coherence( regentlib.exclusive, parts1 ) )
    end
    return coherences
end

return coherence_compute
