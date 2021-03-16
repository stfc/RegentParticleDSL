# Common implementation patterns

This section of the developers guide covers common patters of code used in HartreeParticleDSL, and how they work.

## Data copies between Regions
In many places in HartreeParticleDSL, data is copied from multiple regions, which may be constructured from different
field spaces. Example uses of this functionality is in the HDF5 simple IO module:
```
[mapper_fields:map(function(field)
     return rquote
     particle_array[i].[field.part_field] = hdfreg[i].[field.io_field]
     end
end)];
```
or in the tradequeue implementation in neighbour search:
```
[part_structure:map(function(element)
  return rquote
    tradequeue[int1d(tradequeue_added)].[element.field] = parts[int1d(part)].[element.field]
  end
end)];
```

These implementations make use of Terra's list functionality, some of the utility functions defined in HartreeParticleDSL,
and some other code to achieve automatic generation of data copies between fields.

In this section, we're going to look explicitly at copying from two regions of the same field space, as used in the tradequeue 
implementation. The mapping between regions of different field spaces is similar, but an extra mapping needs to be defined between
the field spaces (as shown in the HDF5 example above, as `field.part_field` and `field.io_field` map to the corresponding fields
in the distinct field spaces).

### Setup
These patterns all start by creating a `terralib.newlist()`, and filling it with a list of fields we want to copy data to and from.
For this example we will use all the fields in our field type, however it could be limited to some predefined subset. We do this by
using the `recurse_fields` and `string_to_fieldpath` utilities:
```
1   local function construct_part_structure()
2     local part_structure = terralib.newlist()
3     local field_strings = {}
4     local type_table = {}
5     for k, v in pairs(part.fields) do
6       recursive_fields.recurse_field(v, field_strings, type_table)
7     end
8     for k, _ in pairs(field_strings) do
9       part_structure:insert({field = string_to_field_path.get_field_path(field_strings[k])})
10    end
11    return part_structure
12  end
13  local part_structure = construct_part_structure()
```
Line 2 creates our Terra list, whilst lines 3 and 4 construct the tables required for the `recursive_fields` utility. We then use the 
`recursive_field` utility on our `part` field space in lines 5 to 7, giving us tables containing all of the (recusively defined) field names
 and types.

In lines 8 to 10, we loop over the field names, and construct the corresponding field paths using the `string_to_field_path` utility. These are
then inserted into the Terra list as a table. In this case the table only contains the field path (accessed through `field`), but we could also choose
to store the type, or a default value, secondary field path etc. as our use case requires.

The Terra list is then returned in line 11, and in line 13 we create a local variable containing the list named `part_structure`.

### Copying between regions.
To implement our copy functionality, we can use the previously constructure Terra list inside a Regent task. Lets imagine we have two identically
sized regions of particles, and want to copy all of the data from `region1` to `region2`, we can create a task:
```
task copy_task(region1 : region(ispace(int1d), part), region2 : region(ispace(int1d), part)) where
  reads(region1), writes(region2) do
  --In this case we're assuming region1 and region2 are identically sized.
  for part in [region1].ispace do
    --Need to write our copy code here.
  end

end
```

We can use the Terra list's `map` function, which creates a map between every element of the list from `A->B`. In this case we want to go from our list (`A`)
to a quote expression containing the code to copy the values between the regions (`B`). This map is called inside an escape as we're generating Regent code using Lua:
```
[part_structure:map( function(element)
  return rquote
    region2[part].[element.field] = region1[part].[element.field]
  end
end)]

```
We define that the `element` variable is each element of the `part_structure` list (the table we input previously). We then create a quote, and access the
two regions at the `part` index, and splice the `.field` field path to access the same field from both regions.

### Design of module imports.

The initial strategy for importing modules was to have everything separate, hence the complete separation of the 2D and 3D periodic tradequeue system. 
However, when designing the non-periodic version, it became apparent that we could create unique 2D and 3D headers, and have each other header required shared between the modules. This is done by the version-specific header setting a global value:
```
DSL_DIMENSIONALITY = 2
```

While this pollutes the global namespace, this can be a reserved variable name. With this set, the shared headers then do:
```
  if DSL_DIMENSIONALITY == 2 then
    neighbour_init = require("src/neighbour_search/cell_pair_tradequeues_nonperiod/2d_neighbour_init")
  elseif DSL_DIMENSIONALITY == 3 then
    neighbour_init = require("src/neighbour_search/cell_pair_tradequeues_nonperiod/3d_neighbour_init")
  end
  return neighbour_init
```
which enables importing only the necessary modules.

If possible, it would be nice to be able to import these headers with a function call, e.g.
```
require("path/to/nonperiodic_header")
set_dimensionality(2)
```
and have that setup the appropriate system, however we have not yet tried implementing this setup.
