# RegentParticleDSL
(Working name)

## Idea
This repository is to cover the goal of building a particle-method DSL using Regent 
metaprogramming. The repository at current contains a few mini examples proving that
this can be done in principle, and will now be fleshed out.

##How does it work
Regent/Terralang/Lua metaprogramming is quite flexible and allows the creation of functions
which essentially inject code into other functions, which can be used to eventually create tasks.
These tasks can be replicated for different injected KERNELS, allowing for a DSL-style approach 
with small sections of user-written code being injected into the larger parallel infrastructure.

##The DSL - defining Kernels
This is still a work in progress. The current idea is that kernels are defined using mostly simple
high level code, though there are a few annoying syntactic hurdles I'd like to avoid that are 
currently required.
A simple pairwise kernel is defined as follows:

    function kernel_one(part1, part2)
    local kernel = rquote
        [part1].extra_variable_1 = 1.0
    end
    return kernel
    end

The function defines `kernel_one` and returns the "injected" code to a higher level function which
will create the task to be executed during the computation. The key part is the section between the
`rquote` and `end`, which controls what code is included in the task.

These kernels can be fed into functions such as:

    local pairwise_task = make_pairwise_interaction(kernel_one, "one")
    local pairwise_task_two = make_pairwise_interaction(kernel_two, "two")

which create the tasks, and the main function can be written by the user to launch a series of tasks
across their particles in whatever order required.


Ideally, I'd like to avoid all the excess "stuff" in the function, and allow definition of inlined 
tasks instead, however at current tasks which have a input of
`task kernel_one(part1 : part, part2 : part)` are unable to actually write to the particle fspace as
Regent is pass by value, so these are copies of some particle, and the only way to do pointers would be
to pass regions, and that ends up too complex for how I envisage this.


##The DSL defining Particle structure
Particles are defined using a field space. There is a `core_part` field space defined in src/particles
which all particles must include, and then the `neighbour_part` field space defined in for 
whichever neighbour search algorith the user wants to use. On top of this, any user-required variables can
be added, as shown in `example_types/example_part.rg`
