#=
Created on 25/01/2023 09:56:39
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

New version of a simulator that uses a graph structure
and should be a bit faster
=#

function compile_circuit(input1, input2, parts...)
    parents = Dict{Int,Int}()
    # find the dependencies of the parts
    for (i, p) in enumerate(parts)
        # check the parent relation
        parents[i] = parent_of(p, parts)
    end
    # find the parents without parents, these are the inputs
    # if acyclic, it is easy to generate the function
    # if there are cycles 


