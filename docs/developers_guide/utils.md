# DEVELOPERS GUIDE: Utility functions (`src/utils`)

Various functionality that powers RegentParticleDSL is designed to be used often throughout the code, and are independent
of the user-level program. These are grouped together as utility functions, and enable various metaprogramming features.

## Compute privileges
The compute privileges file (`src/utils/compute_privileges.rg`) contains functions that analyse a quote (defined by 
`rquote`) and compute the read/write privileges. This is done through analysis of the AST produced by Regent/Terra/Lua
when accesing the table associated with a quote.

### Implementation details - specialized expressions
The `is_specialized_table` object stores a set of ast types that are of interest when doing this AST recursion. From 
testing and looking at `legion/languages/src/regent/ast.t` from the Legion project, we know that all of the required
types are defined by the expressions here, allowing us to make use of the `ast.traverse_node_postorder` functions
defined in regent's `ast` module.

The `is_specialized_node` function takes the table and uses an `ast` function to create a valid way to compare a node
against the `is_specialized_table` object. N.B. This is primarily taken from how Regent uses this function, so 
alternative or improved implementations are welcomed.

### Implementation details - traverse fieldaccess functions.
The `traverse_fieldaccess_postorder_XX_region` functions take an ast node as an input, and search for the sub 
field access elements (`node:is(ast.specialized.expr.FieldAccess)`). If there are sub-FieldAccesses (for example
when accessing nested field spaces) the names of these are appended to give a full set of field names, e.g.
`field1.field2` as well as computing which of the symbols used in the quote they field accesses belong to. There
are implementations of these for one, two or three symbols.

### Implementation details - privilege maps
The `XXX_region_privilege_map` functions take an AST node as an input, and check:
1. If the AST node is an `Assignment` AND the left-hand side of that assignment is a `FieldAccess`, then find the
   field and symbol on which that access took place and add that field to the write table for the relevant symbol.
2. If the AST node is a `FieldAccess`, if its a top-level field and has is applied to one of the symbols, then 
   add it to the read table for the relevant symbol.

IMPORTANT NOTE: Right now for some reason the read accesses are just applied to top-level fields. This could recurse
relatively easily, but doesn't at this time.

### Implementation details - traverse specialized postorder XXX region
The `traverse_specialized_postorder_XXX_region` functions are called on the top level AST node, and for every sub-node
covered by the `is_specialized_node` function, calls the inputted function (a privilege map function).

### Implementation details - XXX region privileges
The `compute_privileges.XXX_region_privileges` functions are called on the kernel function, and return lists of strings, containing the
names of the fields written to and read by the input kernel's quote.

The function generates sufficient symbols to call the kernel, and then creates an instance of the quote, and call the traversal function on
the quote's AST. The function then removes any duplicates from each list, before returning the lists in order:
```
read1, [read2], [read3], write1, [write2], [write3]
```

### Encapsulation
Most of the functions in this utility are not visible outside of the file, with only the `XXX_region_privileges` functions visible externally, 
through the returned table.

## String to fieldpath
The string to fieldpath utility is used in the code to convert field strings (which are easy to abstract from the code) into field path objects, 
which are needed to access fields (particularly recursively defined fields) through metaprogramming in Regent.

### Implementation details - split into table
The `split_into_table` function takes the field string (e.g. `field1.field2`), and splits it into a table of the fields in order from left to right
(so `{field1, field2}`). The `.`s are removed as these are not needed.

### Implementation details get field path
The `string_to_field_path.get_field_path` function takes the field string as an input, and returns the Regent field path associated with that object.
This is done by calling the `split_into_table` function on the input string, then unpacking the resulting table and using it as an input to the 
`regentlib.field_path` function.

### Encapsulation
Only the `string_to_field_path.get_field_path` function is visible externally, through the returned table.

## Recursive fields
The "Recursive fields" utility computes all of the fields and their corresponding types for a field space. For a field space `fspace part` the 
`recursive_fields.recurse_field` function can be used:
```
  local field_strings = {}
  local type_table = {}
  for k, v in pairs(part.fields) do
    recursive_fields.recurse_field(v, field_strings, type_table)
  end
```
and once executed the `field_strings` and `type_table` will contain a full breakdown of the `part` field space. This can be used to set initial
values for all fields, or for copies between field spaces of the same or differing types. A similar concept is used in the HDF5 simple IO module 
to enable the mappers.

### Implementation details
The recurse fields function is fairly straightforward. If the input field's type is a field space instance
(`field.type.is_fspace_instance`) then the function is called recursively on all subfields, passing the name of the field (and any parent fields) 
to the recursive function.

If the input field is not a field space instance, then its name (prepended by any parent field spaces names) is added to the `field_table`, and its
type is added to the `type_table`.
