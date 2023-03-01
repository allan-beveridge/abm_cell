# https://discourse.julialang.org/t/add-new-field-to-struct/73923

using OrderedCollections
using OrdinaryDiffEq
using Distributions
using Random
using CSV
using DataFrames
using Plots
using Printf
#using Setfield

include("parseCell.jl")
#include("Events.jl")
include("PlotCell.jl")

function initialiseEnvironment(envList)
    for env in envList
        ex = Meta.parse(env)
        eval(ex)
    end
    return
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

function convertToFunction(functionName)
    symbol = Symbol(functionName)
    fnc = getfield(Main, symbol)
    return(fnc)
end

function getSymbolList(pairList)
    list = Vector{Symbol}()
    for pair in pairList
        push!(list,pair.first)
    end
    return(list)
end

function convertPairsToFloats(pairList)
    list = []
    for pair in pairList
        push!(list,pair.second)
    end
    return(list)
end

function getCellDictionary(id,position,u0,integrator)
    d = Dict{Symbol,Any}()
    d[Symbol("label__")] = 0
    d[Symbol("id")] = id # + 10
    d[Symbol("pos")] = position
    println("u0:",u0)
    fieldList = Vector{Symbol}()
    for pair in u0
        push!(fieldList,pair.first)
        d[pair.first] = pair.second
    end
    d[Symbol("u0_")] = fieldList
    #d[Symbol("integrator_")] = integrator
    #lineage = Lineage(string(id),1,data.deltat_,"alive",1,Vector{String}())
    #d[Symbol("lineage_")] = lineage
    return(d)
end

Base.@kwdef mutable struct Lineage
    id_::String
    created_::Int32
    deltaT_::Float64
    status_::String # alive, divided, dead
    timeSteps_::Int32
    ancestors_::Vector{String}
end

# @everywhere
mutable struct CellData <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}

    properties::Dict{Symbol,Any}
    lineage_::Lineage
    integrator_
    events_::Vector{Event}
    label_::String
    nSteps_::Int64
end

Base.getproperty(x::CellData, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::CellData, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::CellData) = keys(getfield(CellData, :properties))

# @everywhere
function executeCellEvents(model,cell)
    eventList = getfield(cell,:events_)
    for event in eventList
        name = getfield(event,:name_)
        println("event:",name)
        if(event.test(model,cell,event))
            event.execute(model,cell,event)
        end
    end
end

# @everywhere
function executeModelEvents(model)
    #eventList = getfield(model,:events)
    eventList = model.events
    cell = nothing
    for event in eventList
        name = getfield(event,:name_)
        println("executeModelEvents:",name)
        if(event.test(model,cell,event))
            event.execute(model,cell,event)
        end
    end
end

##@everywhere 
function update_cell!(agent, model)
    println("update_cell!")
    #println("cell:",agent)
    nSteps = getfield(agent,:nSteps_) + 1
    setfield!(agent,:nSteps_,nSteps)
    integrator = getfield(agent,:integrator_)
    OrdinaryDiffEq.step!(integrator, model.dt, true)
    index = 1
    for field in agent.u0_
        setproperty!(agent,field,integrator.u[index])
        index += 1
    end
    executeCellEvents(model,agent)
    return()
end

##@everywhere 
function update_model!(model) #2
    model.nSteps += 1
    println("update_model!")
    println(model)
    executeModelEvents(model)
    return
    updateNutrients(model)
    return
end

# @everywhere
function generateCellModel(data)
    prob = ODEProblem(data.rn_, data.u0_, data.tspan_, data.p_)
    integrator = OrdinaryDiffEq.init(prob, OrdinaryDiffEq.Rodas4(); advance_to_tstop = true)

    model = ABM(
        CellData; #Cell1;
        properties = Dict(
            :rn => data.rn_,
            :p => data.p_,
            :u0 => data.u0_,
            :nutrients => data.nutrients_,
            :events => nothing,
            :tspan => data.tspan_,
            :lastCell => data.nCells_,
            :nAgents => data.nCells_,
            :i => integrator, # The OrdinaryDiffEq integrator,
            :dt => data.deltat_,
            :nSteps => 0
        ),
    )
    model.events = getModelEvents(data.nutrients_)
    println("nutrients:",data.nutrients_)
    u0 = convertPairsToFloats(data.u0_)
    println(u0)
    println(typeof(u0))
    println(data.u0_)
    println(typeof(data.u0_))

    fieldList = getSymbolList(data.u0_)
    println(fieldList)
    println(typeof(fieldList))
    uList = nothing
    events = getCellEvents()
    cellLabel = "label"
    for i in 1:model.nAgents
        position = (i*10.0,0.0)
        prob = ODEProblem(model.rn, u0, data.tspan_, data.p_)
        integrator = OrdinaryDiffEq.init(prob, OrdinaryDiffEq.Rodas4(); advance_to_tstop = true)
        lineage = Lineage(string(i),1,data.deltat_,"alive",1,Vector{String}())
        d = getCellDictionary(i,position,data.u0_,integrator)
        #d[:label__] = 2
        d[Symbol("label__")] = 2
# @everywhere
        cell = add_agent!(position,model,d,lineage,integrator,events,cellLabel,0)
        cell.LABEL = "LAbel"
        setproperty!(cell,:label__,3)
        println("cell:",cell)
        println("LABEL:",cell.LABEL)
        println(cell.a)
    end
    return(model)
end


