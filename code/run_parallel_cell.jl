using Dates
using Distributed
@everywhere using Agents
#using Random

include("CommandParser.jl")
@everywhere begin
    include("Cell.jl")
end

# https://discourse.julialang.org/t/using-local-module-with-everywhere/87163
# julia -p 2 run_parallel_cell.jl -cpus 2

@everywhere function randomSeed()
    date = now()
    cols = split(string(date),":")
    time = replace(cols[length(cols)],"." => "")
    time = parse(Int32,time)
    println("random:",time)
    Random.seed!(time)
    return
end

function runParallelCell(modelParams,outDir,cpus)
    println("cpus:",cpus, " np:",nprocs(), " nt:", Base.Threads.nthreads())
    if(nprocs() == 1)
        println("only 1 cpu, ...aborting")
        exit()
    end

    maxTime = modelParams.tspan_[2]
    nCycles = Int64(modelParams.tspan_[2]/modelParams.deltat_)
    nCycles = 200
    modelList = Vector{ABM}()
    for i in 1:cpus
        outputDir = outDir * "_" * string(i)
        println("o:",outputDir)
        randomSeed() 
        model = generateCellModel(modelParams,i*1000)
        push!(modelList,model)
    end
    println("runParallelCell:",cpus)
    aList = [:r,:e_t,:e_m,:q,:m_r,:m_t,:m_m,:m_q,:c_r,:c_t,:c_m,:c_q,:a,:s_i]
    aList = [:id_,:r,:e_t,:e_m,:q,:m_r,:m_t,:m_m,:m_q,:c_r,:c_t,:c_m,:c_q,:a,:s_i]
    agentResults, modelResults = ensemblerun!(modelList, update_cell!, update_model!, nCycles; mdata = [], adata = aList, parallel=true)

    CSV.write("p_abm_cell.csv",agentResults)

    roundColumns(agentResults,aList)
    CSV.write("parallel_abm_cell.csv",agentResults)
    plotPrefix = "agent_model_dt=" * string(modelParams.deltat_) * "_maxt=" * string(maxTime)
    println("prefix:",plotPrefix)
    model = modelList[1]
    agentPlot(agentResults,model.nAgents,outDir,plotPrefix)
    exit()

    eventList = model.events
    for event in eventList
        println("event:",event)
        event.save(event)
    end
    return
end

function run(program,dataFile,outputDir,cpus)
    createDirectory(outputDir)
    modelParams = parseCell(dataFile) 
    initialiseEnvironment(modelParams.envList_) 
    println("env:",ENV)
    println("cpus:",cpus, " ", typeof(cpus))
    if(cpus > 1)
        runParallelCell(modelParams,outputDir,cpus)
    end
    return()
end

function main()
    dataFile = "m.dat"
    outputDir = "./parallel"
    parameterDictionary = CommandParser.readParameters("run_cell.dat")
    println(parameterDictionary)
    dict = CommandParser.parse(parameterDictionary,ARGS)
    #println(dict)

    program = nothing
    cpus = 1
    if (haskey(dict,"-run"))
       program = dict["-run"]
    end
    if (haskey(dict,"-cpus"))
       cpus = dict["-cpus"]
    end
    run(program,dataFile,outputDir,cpus)
    return
end

main()

