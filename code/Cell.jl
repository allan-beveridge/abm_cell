# https://discourse.julialang.org/t/add-new-field-to-struct/73923

using OrderedCollections
#using Agents
using OrdinaryDiffEq
using Distributions
using Random
using CSV
using DataFrames
using Plots
using Printf
#using Setfield

include("parseCell.jl")
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

function getCellDictionary(agent_id,id_,position,u0,integrator)
    #println("getCellDictionary:",id_)
    d = Dict{Symbol,Any}()
    d[Symbol("label__")] = 0
    d[Symbol("id")] = agent_id
    d[Symbol("id_")] = id_ 
    d[Symbol("pos_")] = position
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

function executeCellEvents(model,cell)
    #println("executeCellEvents")
    eventList = getfield(cell,:events_)
    for event in eventList
        name = getfield(event,:name_)
        if(event.test(model,cell,event))
            event.execute(model,cell,event)
        end
    end
    return
end

function executeModelEvents(model)
    if(model.events == nothing)
        return
    end
    eventList = model.events
    cell = nothing
    for event in eventList
        name = getfield(event,:name_)
        #println("executeModelEvents:",name)
        if(event.test(model,cell,event))
            event.execute(model,cell,event)
        end
    end
    return
end

function no_events_update_cell!(agent, model)
    nSteps = getfield(agent,:nSteps_) + 1
    setfield!(agent,:nSteps_,nSteps)
    integrator = getfield(agent,:integrator_)
    OrdinaryDiffEq.step!(integrator, model.dt, true)
    index = 1
    for field in agent.u0_
        setproperty!(agent,field,integrator.u[index])
        index += 1
    end
    return
end

function no_events_update_model!(model) #2
    model.nSteps += 1
    return
end

function update_cell!(agent, model)
    #println("update_cell!")
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
    #println("END update_cell!")
    return()
end

function update_model!(model) #2
    model.nSteps += 1
    #println("Update_model!")
    executeModelEvents(model)
    #updateNutrients(model)
    return
end

function set_agent_id(cell,id)
    setfield!(cell,:id,id)
end

function set_cell_id(cell,id)
    setproperty!(cell,:id,id)
end

function generateCellModel(data,cellIndex=0)
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
        cell_id = cellIndex + i
        d = getCellDictionary(i,cell_id,position,data.u0_,integrator)
        #d[:label__] = 2
        d[Symbol("label__")] = 2
        cell = add_agent!(position,model,d,lineage,integrator,events,cellLabel,0)
        cell.rn_ = data.rn_
        cell.p_ = data.p_
        cell.LABEL = "LAbel"
        setproperty!(cell,:label__,3)
        #set_agent_id(cell,666)
        #set_cell_id(cell,555)
        println("LABEL:",cell.LABEL)
        println(cell.a)
        #dump(cell)
    end
    return(model)
end

function xxx_runCell(dataFile,outputDir)
    ##ENV["GKSwstype"] = "100"
    Random.seed!(6549)
    createDirectory(outputDir)
    model_params = parseCell(dataFile)
    initialiseEnvironment(model_params.envList_)
    #println("env:",ENV)
    model = generateCellModel(model_params)

    maxTime = model_params.tspan_[2]
    nCycles = Int64(maxTime/model_params.deltat_)
    nCycles = 10
    nCycles = 500
    nCycles = 200
    aList = [:r,:e_t,:e_m,:q,:m_r,:m_t,:m_m,:m_q,:c_r,:c_t,:c_m,:c_q,:a,:s_i]
    aList = [:id_,:r,:e_t,:e_m,:q,:m_r,:m_t,:m_m,:m_q,:c_r,:c_t,:c_m,:c_q,:a,:s_i]
    agentResults, modelResults = run!(model, update_cell!, update_model!, nCycles; mdata = [], adata = aList)

    roundColumns(agentResults,aList)
    CSV.write("abm_cell.csv",agentResults)
    plotPrefix = "agent_model_dt=" * string(model_params.deltat_) * "_maxt=" * string(maxTime)
    println("prefix:",plotPrefix)
    agentPlot(agentResults,model.nAgents,outputDir,plotPrefix)

    eventList = model.events
    for event in eventList
        println("event:",event)
        event.save(event)
    end

    println("finished")

end


#dataFile = "m.dat"
#outputDir = "./model"
#runCell(dataFile,outputDir)

