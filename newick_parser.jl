import Automa
import Automa.RegExp: @re_str
const re = Automa.RegExp

machine = (function ()
    # Primitives
    label      = re"[^\(\);,:]+"
    len        = re":[^\(\);,:]+"
    cladestart = re"\("
    cladestop  = re"\)"
    sep        = re","
    finish     = re";"

    # Composites
    metadata   = (label * len) | label | len
    nodechange = cladestart | cladestop | sep
    element    = nodechange * re.opt(metadata)
    tree       = re.opt(metadata) * re.rep(element) * finish

    # Actions
    label.actions[:enter] = [:mark]
    label.actions[:exit]  = [:label]
    len.actions[:enter] = [:mark]
    len.actions[:exit] = [:len]
    cladestart.actions[:enter] = [:push]
    cladestop.actions[:enter] = [:pop]
    sep.actions[:enter] = [:pop, :push]
    finish.actions[:enter] = [:finish]

    # Finally, compile the final FASTA pattern into a state machine.
    return Automa.compile(tree)
end)()

mutable struct Node
    label::Union{Nothing, String}
    len::Union{Nothing, Float64}
    parent::Union{Nothing, Node}
end

function Base.show(io::IO, node::Node)
    parent = node.parent === nothing ? nothing : "Node(\"$(node.parent.label)\")"
    print(io, "Node(\"$(node.label)\", $(node.len), $parent)")
end

@noinline function popstack!(stack::Vector{Node})
    length(stack) > 1 || error("Too many closing parenthesis in file")
    pop!(stack)
    return nothing
end

actions = Dict(
    :finish => :(finished = true),
    :mark => :(mark = p),
    :label => :(last(stack).label = String(data[mark:p-1])),
    :len => :(last(stack).len = Base.parse(Float64, String(data[mark+1:p-1]))),
    :pop => :(popstack!(stack)),
    :push => quote
        push!(stack, Node(nothing, nothing, last(stack)))
        push!(nodes, last(stack))
    end
)

context = Automa.CodeGenContext()
@eval function parse(data::Union{String,Vector{UInt8}})
    mark = 0
    nodes = [Node(nothing, nothing, nothing)]
    stack = copy(nodes)
    finished = false

    $(Automa.generate_init_code(context, machine))
    p_end = p_eof = lastindex(data)
    $(Automa.generate_exec_code(context, machine, actions))

    if (cs != 0) & !finished
        error("failed to parse on byte ", p)
    end

    # Final stack should only have the root node.
    if length(stack) != 1
        error("Not enough closing parentheses in file")
    end

    return nodes
end
