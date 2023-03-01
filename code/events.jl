function alwaysUpdate(model,cell,event)
    return(true)
end

function dontSave(model,cell,event)
    return()
end

#Base.getproperty(x::Event, property::Symbol) = getfield(x, :functions_)[property]
#Base.setproperty!(x::Event, property::Symbol, value) = getfield(x, :functions_)[property] = value
#Base.propertynames(x::Event) = keys(getfield(Event, :functions_))

# new event code

function move(model,agent,event)
    fsym, compression, torque = transform_forces(agent)
    direction =  model.dt * model.mobility .* fsym
    walk!(agent, direction, model)
    agent.length += model.dt * model.mobility .* compression
    agent.orientation += model.dt * model.mobility .* torque
    agent.growthprog += model.dt * agent.growthrate
    update_nodes!(agent)
    return agent.pos
end

function divideByGrowth(model,agent,event)
    if a.growthprog ≥ 1
        return(true)
    end
    return(false)
end

function divideCells(model,agent,event)
    model.stepCount += 1
    for a in allagents(model)
        if a.growthprog ≥ 1
            ## When a cell has matured, it divides into two daughter cells on the
            ## positions of its nodes.

            id = model.lastCell + 1
            d = getCellDictionary(id,a.p1,0.0,a.orientation,0.0,0.1 *rand(model.rng) + 0.05)
            eventList = getfield(a,:events_)
            events::Vector{Event} = deepcopy(eventList)
            cell = add_agent!(a.p1,model,d,events)
            update_nodes!(cell)

            id = model.lastCell + 2
            d = getCellDictionary(id,a.p2,0.0,a.orientation,0.0,0.1 *rand(model.rng) + 0.05)
            events = deepcopy(eventList)
            cell = add_agent!(a.p2,model,d,events)
            update_nodes!(cell)
            model.lastCell += 2

            kill_agent!(a, model)
        else
            ## The rest lengh of the internal spring grows with time. This causes
            ## the nodes to physically separate.
            uv = unitvector(a.orientation)
            internalforce = model.hardness * (a.length - a.growthprog) .* uv
            a.f1 = -1 .* internalforce
            a.f2 = internalforce
        end
    end
    ## Bacteria can interact with more than on other cell at the same time, therefore,
    ## we need to specify the option `:all` in `interacting_pairs`
    for (a1, a2) in interacting_pairs(model, 2.0, :all)
        interact!(a1, a2, model)
    end
    return
end

#########################################################################

# Some geometry convenience functions
unitvector(φ) = reverse(sincos(φ))
cross2D(a, b) = a[1] * b[2] - a[2] * b[1]
#nothing # hide

function update_nodes!(a)
    println("update_nodes")
    offset = 0.5 * a.length .* unitvector(a.orientation)
    a.p1 = a.pos .+ offset
    a.p2 = a.pos .- offset
end

function transform_forces(agent)
    #println("transform_forces")
    ## symmetric forces (CM movement)
    fsym = agent.f1 .+ agent.f2
    ## antisymmetric forces (compression, torque)
    fasym = agent.f1 .- agent.f2
    uv = unitvector(agent.orientation)
    compression = dot(uv, fasym)
    torque = 0.5 * cross2D(uv, fasym)
    return fsym, compression, torque
end

function update_nodes!(a)
    println("update_nodes")
    offset = 0.5 * a.length .* unitvector(a.orientation)
    a.p1 = a.pos .+ offset
    a.p2 = a.pos .- offset
end

# ### Helper functions
function interact!(a1, a2, model)
    #println("interact!")
    n11 = noderepulsion(a1.p1, a2.p1, model)
    n12 = noderepulsion(a1.p1, a2.p2, model)
    n21 = noderepulsion(a1.p2, a2.p1, model)
    n22 = noderepulsion(a1.p2, a2.p2, model)
    a1.f1 = @. a1.f1 + (n11 + n12)
    a1.f2 = @. a1.f2 + (n21 + n22)
    a2.f1 = @. a2.f1 - (n11 + n21)
    a2.f2 = @. a2.f2 - (n12 + n22)
end

function noderepulsion(p1::NTuple{2,Float64}, p2::NTuple{2,Float64}, model::ABM)
    #println("noderepulsion")
    delta = p1 .- p2
    distance = norm(delta)
    if distance ≤ 1
        uv = delta ./ distance
        return (model.hardness * (1 - distance)) .* uv
    end
    return (0, 0)
end


