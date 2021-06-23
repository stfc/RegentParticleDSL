# Writing Efficient Kernels

While HartreeParticleDSL aims to make things as easy as possible to get performant particle simulations, there are a few small but important ways to improve
performance.

## Avoid Symmetric Pairwise kernels
At the moment, the Symmetric Pairwise kernels are not yet supported by the high performance sections of the runtime system. We plan to add this when we can, 
however if this particularly important for your use-case, please comment on the issue at https://github.com/stfc/RegentParticleDSL/issues/92 .

## Use reduction operators

When writing kernels for particle pairs, parallelism will be restricted if a particle is updated using the `=` operator. To avoid this issue, particle updates
should be done with "reduction" operators wherever possible, i.e. `+=`, `-=`, `*=`, `/=`, `max=` or `min=`. 

To use max or min operators, there are special inbuilt functions for `max` and `min`, e.g.
```
part1.max_visc max= max(option1, option2)
```

If the code is instead `max(part1.max_visc, visc)`, this can be written as:
```
part1.max_visc max= visc
```

## Avoiding reading and writing/reduction from the same fields in a single kernel
When writing kernels, the code will be able to do a better job of optimising performance if fields are read-only or write/reduce-only. 
This restriction is lessened for per-particle tasks, reading and writing from the same field should not impact performance in that case.
