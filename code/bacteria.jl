# https://juliadynamics.github.io/Agents.jl/v3.3/examples/growing_bacteria/

using Agents, LinearAlgebra
using Random
using CSV

using InteractiveDynamics
using CairoMakie # choose plotting backend

include("events.jl")
include("inputParser.jl")

mutable struct zzz_Event
    name_::String
    data_::Dict{String,Any}
    functions_::Dict{Symbol,Any}
    global_::Bool
end

Base.getproperty(x::Event, property::Symbol) = getfield(x, :functions_)[property]
Base.setproperty!(x::Event, property::Symbol, value) = getfield(x, :functions_)[property] = value
Base.propertynames(x::Event) = keys(getfield(Event, :functions_))

mutable struct Cell <: AbstractAgent
    id::Int
    pos::Tuple{Float64, Float64}

    properties::Dict{Symbol,Any}
    events_::Vector{Event}
end

Base.getproperty(x::Cell, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::Cell, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::Cell) = keys(getfield(Cell, :properties))

function cassini_oval(agent)
    t = LinRange(0, 2π, 50)
    a = agent.growthprog
    b = 1
    m = @. 2 * sqrt((b^4 - a^4) + a^4 * cos(2 * t)^2) + 2 * a^2 * cos(2 * t)
    C = sqrt.(m / 2)

    x = C .* cos.(t)
    y = C .* sin.(t)

    uv = reverse(sincos(agent.orientation))
    θ = atan(uv[2], uv[1])
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]

    bacteria = R * permutedims([x y])
    coords = [Point2f(x, y) for (x, y) in zip(bacteria[1, :], bacteria[2, :])]
    scale(Polygon(coords), 0.5)
end

######################################################################

function getCellEvents()::Vector{Event}
    eventList = Vector{Event}()
    data = Dict{String,Any}()
    functions = Dict{Symbol,Any}()
    functions[:test] = alwaysUpdate
    functions[:execute] = move
    functions[:save] = dontSave
    event = Event("move",data,functions)
    push!(eventList,event)
    return(eventList)
end

function getModelEvents()
    eventList = Vector{Event}()
    data = Dict{String,Any}()
    functions = Dict{Symbol,Any}()
    functions[:test] = alwaysUpdate
    functions[:execute] = divideCells
    functions[:save] = dontSave
    event = Event("divide",data,functions)
    push!(eventList,event)
    return(eventList)
end

function update_cell!(agent::Cell, model)
    println("update_cell")
    event = nothing
    eventList = getfield(agent,:events_)
    for event in eventList
        name = getfield(event,:name_)
        if(event.test(model,agent,event))
            event.execute(model,agent,event)
        end
    end
    return
end

function update_model!(model)
    println("update_model:",model.stepCount)
    model.stepCount += 1
    for event in model.events_
        event.execute(model,nothing,event)
    end
    return
end

function getCellDictionary(varMap,varTypeMap)
    d = Dict{Symbol,Any}()
    for (key,value) in varMap
        d[Symbol(key)] = value
    end

    p1::NTuple{2,Float64} = (0.0,0.0)
    p2::NTuple{2,Float64} = (0.0,0.0)
    f1::NTuple{2,Float64} = (0.0,0.0)
    f2::NTuple{2,Float64} = (0.0,0.0)

    orientation = d[:orientation]
    position = d[:pos]
    length = d[:length]
    offset = 0.5 * length .* unitvector(orientation)
    p1 = position[1] .+ offset
    p2 = position[2] .- offset

    d[Symbol("p1")] = p1
    d[Symbol("p2")] = p2
    d[Symbol("f1")] = f1
    d[Symbol("f2")] = f2

    return(d)
end

#function getCellDictionary(agent_id,position,length,orientation,growthprog,growthrate)
function getCellDictionary(agent_id::Int,position::NTuple{2,Float64},length::Float64,orientation::Float64,growthprog::Float64,growthrate::Float64)

    pos::NTuple{2,Float64} = (position[1],position[2])
    d = Dict{Symbol,Any}()
    d[Symbol("id")] = agent_id
    d[Symbol("pos")] = pos

    d[Symbol("length")] = length
    d[Symbol("orientation")] = orientation
    d[Symbol("growthprog")] = growthprog
    d[Symbol("growthrate")] = growthrate

    p1::NTuple{2,Float64} = (0.0,0.0)
    p2::NTuple{2,Float64} = (0.0,0.0)
    f1::NTuple{2,Float64} = (0.0,0.0)
    f2::NTuple{2,Float64} = (0.0,0.0)

    #println("update_nodes")
    #offset = 0.5 * length .* unitvector(orientation)
    #p1 = position[1] .+ offset
    #p2 = position[2] .- offset

    d[Symbol("p1")] = p1
    d[Symbol("p2")] = p2
    d[Symbol("f1")] = f1
    d[Symbol("f2")] = f2

    return(d)
end

function createDirectory(directory)
    if !isdir(directory)
        path = mkdir(directory)
    else
       println("Warning, directory already exists:" * directory)
    end
    return
end

function roundColumns(df,columnList,nDigits=3)
    for column in columnList
        name = Symbol(column)
        df[:,name] = round.(df[:,name];digits=nDigits)
    end
end

function roundTuples(tupleList,nDigits=3)
    list = []
    for tuple in tupleList
        t = (round.(tuple[1];digits=nDigits),round(tuple[2];digits=nDigits))
        push!(list,t)
    end
    return(list)
end

function roundDataFrame(df,columnList,nDigits=3)
    for col in columnList
        column = df[:,col]
        tupleList = roundTuples(column,nDigits)
        df[:,col] = tupleList
    end
    return
end

function getEvents(eventsMap,eventDataList)

    println("getEvents")
    matchedEventsMap = OrderedDict{String,Event}()
    eventList = []
    for data in eventDataList
        eventLabel = data[1]
        eventName = data[2]
        global_ = data[3]
        variableList = data[4]
        if(!haskey(eventsMap,eventName))
            abort__("Error, invalid event:-1" * eventName)
        end
        if(haskey(matchedEventsMap,eventName))
            currentEvent = matchedEventsMap[eventName]
        else
            currentEvent = eventsMap[eventName]
        end
        if(global_)
            event = currentEvent
        else
            event = deepcopy(currentEvent)
        end
        resetEvent(event,variableList)
        matchedEventsMap[eventLabel] = event
    end
    return(matchedEventsMap)
end

function old_getEvents(label,eventsMap,cellEventsMap)
    events = cellEventsMap[label]
    for event_name in events
        println(event_name)
        if( !haskey(eventsMap,event_name) )
            abort__("Error, invalid event:" * event_name)
        end
        event = eventsMap[event_name]
        push!(eventList,event)
    end
    return(eventList)
end

function parseCells(model,variableMap,cellList,eventList)
    for line in cellList
        cols = split(strip(line))
        id = strip(cols[3])
        x = Base.parse(Float64,strip(cols[4]))
        y = Base.parse(Float64,strip(cols[5]))
        length = Base.parse(Float64,strip(cols[6]))
        orientation = Base.parse(Float64,strip(cols[7]))
        growthprog = Base.parse(Float64,strip(cols[8]))
        growthrate = Base.parse(Float64,strip(cols[9]))
        position = (x,y)
        map = Dict{String,Any}("id" => 1, "pos" => position, "length" => length, "orientation" => orientation, "growthprog" => growthprog, "growthrate" => growthrate)
        d = getCellDictionary(map,variableMap)
        cellEventList = deepcopy(eventList)
        cell = add_agent!(position,model,d,cellEventList)
        println(d)
    end
    return
end

function initialiseCells(model,variableMap,eventList)
    position = (6.5,4.0)
    map = Dict{String,Any}("id" => 1, "pos" => position, "length" => 0.0, "orientation" => 0.3, "growthprog" => 0.0, "growthrate" => 0.1)
    d = getCellDictionary(map,variableMap)
    cellEventList = deepcopy(eventList)
    cell = add_agent!(position,model,d,cellEventList)

    position = (7.5,4.0)
    map = Dict{String,Any}("id" => 2, "pos" => position, "length" => 0.0, "orientation" => 0.0, "growthprog" => 0.0, "growthrate" => 0.1)
    d = getCellDictionary(map,variableMap)
    cellEventList = deepcopy(eventList)
    eventList = deepcopy(cellEventList)
    cell = add_agent!(position,model,d,eventList)
    return
end

function demo(inputFile)
    println("reading:",inputFile)
    dataMap = parseInputFile(inputFile)
    modelEventsDataList = dataMap["modelEventsList"]
    eventsMap = deepcopy(dataMap["EventsMap"])
    cellEventsMap = deepcopy(dataMap["cellEventsMap"])
    cellEventDataList = cellEventsMap["SimpleCell"]
    map = getEvents(eventsMap,cellEventDataList) 
    cellEventList = collect(values(map))

    modelEventsMap = deepcopy(dataMap["modelEventsMap"])
    map = getEvents(modelEventsMap,modelEventsDataList) 
    modelEventList = collect(values(map))
println(modelEventList)
    
    space = ContinuousSpace((14, 9), 1.0; periodic = false)
    model = ABM(
        Cell,
        space,
        properties = Dict(:dt => 0.005, :hardness => 1e2, :mobility => 1.0, :stepCount => 0, :lastCell => 2, :events_ => Vector{Event}()),
        rng = MersenneTwister(1680)
    )
    model.events_ = modelEventList

    variableMap = dataMap["SimpleCell"]

    #parseCells(model,variableMap,dataMap["cellList"],cellEventList)
    initialiseCells(model,variableMap,cellEventList)

    adata = [:pos, :length, :orientation, :growthprog, :p1, :p2, :f1, :f2]

    CairoMakie.activate!() # hide

##### run #####
    nCycles = 5000
    nCycles = 200
    nCycles = 2100
    nCycles = 2001
    nCycles = 20
    nCycles = 2002

    println("run")
    aList = [:pos,:orientation,:length,:growthprog,:growthrate,:p1,:p2,:f1,:f2]
    cellResults, modelResults = run!(model, update_cell!, update_model!, nCycles; mdata = [], adata = aList)

    outputDir = "bacteria"
    csvFile = outputDir * "/bacteria.csv"
    createDirectory(outputDir)

    aList = [:orientation,:length,:growthprog,:growthrate]
    roundColumns(cellResults,aList,3)
    columnList = ["pos","p1","p2","f1","f2"]
    roundDataFrame(cellResults,columnList)

    println(csvFile)
    CSV.write(csvFile,cellResults)
    println("finished")
    return
end

function main()
    inputFile = "bacteria.dat"
    demo(inputFile)
    return
end

main()
