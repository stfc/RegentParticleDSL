# Example templates

The `templates` directory contains templates that can be used to as a basis for creating new code to use with the
DSL. Each file contains the basis of a separate part of the DSL.

## Particle templates
The `templates/part_templates.rg` file contains all the essential boilerplate code required to create a new particle
type. This file has two imports. The first import:
```
require("src/particles/core_part")
```
imports the `core_part` type, and is required for all particle types.

The second import:
```
require("src/neighbour_search/cell_pair/cell_pair")
```
is specific to the neighbour search algorithm chosen, in this example the periodic cell pair algorithm. For
information on what import to use for a specific neighbour search algorithm, check the appropriate documentation.

The rest of the file details the particle type itself. The declaration of the particle type and the `core_part_space` and
`neighbour_part_space` must exactly follow the template, as these are visible througout the DSL.

The `extra_variable` fields are optional, and can be named as required. If your particle method has any data not
covered by the `core_part` type, it can be added here, and is visible to any kernels you create.

For a full list of the fields contained in the `core_part` check the [main DSL structure documentation](DSL.md), or 
view the declaration in `src/particles/core_part.rg`.

## Kernel templates
The `templates/kernel_templates.rg` contains templates for the three main types of Kernel's currently supported in 
RegentParticleDSL:
1. Pairwise particle kernels.
2. Asymmetric pairwise particle kernels.
3. Per-particle kernels.

For all of these templates, almost all of the boilerplate code (names/variable names etc.) is adjustable, other than the 
Regent/Lua keywords (`function`, `local`, `rquote`, `end`, `return`).
