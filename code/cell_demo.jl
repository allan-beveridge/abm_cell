include("inputParser.jl")
include("demo2_events.jl")
include("CommandParser.jl")

function plotCellCounts(cellCountFile,plotFile)
    df = DataFrame(CSV.File(cellCountFile))
    plot(df.timeStep,df.cellCount,xlabel="Time",ylabel="Cell Count",color="blue")
    savefig(plotFile)
    return
end

function initialiseModel__(abm_data,cellIndex=0)

    println("initialiseModel__")
    cellMap = abm_data["cellMap"]
    data = cellMap["A"]

    prob = ODEProblem(data.rn_, data.u0_, data.tspan_, data.p_)
    integrator = OrdinaryDiffEq.init(prob, OrdinaryDiffEq.Rodas4(); advance_to_tstop = true)

    #tspan__ = (0,4000)

    model = ABM(
        CellData; #Cell1;
        properties = Dict(
            :rn => data.rn_,
            :p => data.p_,
            :u0 => data.u0_,
            :cellDataMap => cellMap,
            :events => nothing,
            :tspan => data.tspan_,
            :lastCell => data.nCells_,
            :nAgents => data.nCells_,
            :i => integrator, # The OrdinaryDiffEq integrator,
            :dt => data.deltat_,
            :nSteps => 0
        ),
    )
    #u0 = convertPairsToFloats(data.u0_)
    #fieldList = getSymbolList(data.u0_)

    initialiseCells__(model,abm_data,cellIndex)

    modelEventsMap =  getModelEvents(abm_data)
    model.events =  collect(values(modelEventsMap))
    return(model)

end

function initialiseCells__(model,abm_data,cellIndex=0)
    cellMap = abm_data["cellMap"]
    cellList = parseCells(abm_data["cellList"])
    u0_map = abm_data["u0"]
    nCells = length(cellList)
    cell_u0 = nothing
    #for i in 1:model.nAgents
    for i in 1:nCells
        cellType, cell_id, position = cellList[i]
        println(cellType, " ", cell_id, " ", position)
        cellData = cellMap[cellType]

        events = getCellEvents(abm_data,cellType)

        if(haskey(u0_map,cellType))
            cell_u0 = deepcopy(u0_map[cellType])
        else
            cell_u0 = deepcopy(cellData.u0_)
        end
        u0 = convertPairsToFloats(cell_u0)

        modelParams = cellMap[cellType]
        prob = generate_ode_problem(modelParams,u0,cellData.tspan_)
        #prob = ODEProblem(model.rn, u0, cellData.tspan_, cellData.p_)
        integrator = OrdinaryDiffEq.init(prob, OrdinaryDiffEq.Rodas4(); advance_to_tstop = true)
        lineage = Lineage(string(i),1,cellData.deltat_,"alive",1,Vector{String}())
        cell_id = cellIndex + i
        d = getCellDictionary(i,cell_id,position,cell_u0,integrator)
        #d = getCellDictionary(i,cell_id,position,cellData.u0_,integrator)
        d[Symbol("label__")] = 2
        d[Symbol("label")] = 2
        cell = add_agent!(position,model,d,lineage,integrator,events,cellType,0)
        cell.p_ = cellData.p_
        cell.rn_ = cellData.rn_
    end
    model.nAgents = nCells
    model.lastCell = nCells
    return()
end

function generate_ode_problem(modelParams,u0,tspan)
    prob = ODEProblem(modelParams.rn_, u0, tspan, modelParams.p_)
    return(prob)
end


function testIntegrate__(abm_data,outPrefix,cellIndex=0)
    Random.seed!(6549)
    outputDir, outputPrefix = abm_data["output"]
    createDirectory(outputDir)
    outputDir = abspath(outputDir)
    println("O:",outputDir, " ", outputPrefix)
    model = initialiseModel__(abm_data,cellIndex)

    cellMap = abm_data["cellMap"]
    data_a = cellMap["A"]
    modelParams = data_a
    maxTime = modelParams.tspan_[2]
    nCycles = Int64(maxTime/modelParams.deltat_)
    nCycles = abm_data["nCycles"]
    #nCycles = 400000
    println("nCycles:",nCycles)
    run_info = "dt=" * string(modelParams.deltat_) * "_maxt=" * string(maxTime)
    aList = [:r,:e_t,:e_m,:q,:m_r,:m_t,:m_m,:m_q,:c_r,:c_t,:c_m,:c_q,:a,:s_i]
    aList = [:id_,:r,:e_t,:e_m,:q,:m_r,:m_t,:m_m,:m_q,:c_r,:c_t,:c_m,:c_q,:a,:s_i]
    agentResults, modelResults = run!(model, no_events_update_cell!, no_events_update_model!, nCycles; mdata = [], adata = aList)

    roundColumns(agentResults,aList)
    csvFile = outputDir * "/" * outPrefix * "_general_abm_cell.csv"
    CSV.write(csvFile,agentResults)
    println(csvFile)

    plotPrefix = outPrefix * "_" * run_info
    println("prefix::",plotPrefix)
    agentPlot(agentResults,model.nAgents,outputDir,plotPrefix)
    return

    saveCellEvents(model)
    eventList = model.events
    for event in eventList
        event.save(event)
    end
    println("finished")
    return
end

function run__(abm_data,outPrefix,cellIndex=0)
    command = "export GKSwstype=\"100\""
    ex = Meta.parse(command)
    #eval(ex)
    Random.seed!(6549)
    outputDir, outputPrefix = abm_data["output"]
    createDirectory(outputDir)
    outputDir = abspath(outputDir)
    println("O:",outputDir, " ", outputPrefix)
    model = initialiseModel__(abm_data,cellIndex)

    cellMap = abm_data["cellMap"]
    data_a = cellMap["A"]
    modelParams = data_a
    maxTime = modelParams.tspan_[2]
    nCycles = Int64(maxTime/modelParams.deltat_)
    nCycles = abm_data["nCycles"]
    nCycles = 1000
    println("nCycles:",nCycles)
    run_info = "dt=" * string(modelParams.deltat_) * "_maxt=" * string(maxTime)
    aList = [:r,:e_t,:e_m,:q,:m_r,:m_t,:m_m,:m_q,:c_r,:c_t,:c_m,:c_q,:a,:s_i]
    aList = [:id_,:r,:e_t,:e_m,:q,:m_r,:m_t,:m_m,:m_q,:c_r,:c_t,:c_m,:c_q,:a,:s_i]
    #aList = [:label__,:r,:e_t,:e_m,:q,:m_r,:m_t,:m_m,:m_q,:c_r,:c_t,:c_m,:c_q,:a,:s_i]
    agentResults, modelResults = run!(model, update_cell!, update_model!, nCycles; mdata = [], adata = aList)
    #CSV.write("zzz_abm_cell.csv",agentResults)

    roundColumns(agentResults,aList)
    csvFile = outputDir * "/" * outPrefix * "_general_abm_cell.csv"
    CSV.write(csvFile,agentResults)
    println(csvFile)

    #plotPrefix = "agent_model_dt=" * string(modelParams.deltat_) * "_maxt=" * string(maxTime)
    plotPrefix = outPrefix * "_" * run_info
    println("prefix::",plotPrefix)
    agentPlot(agentResults,model.nAgents,outputDir,plotPrefix)

    saveCellEvents(model)
    eventList = model.events
    for event in eventList
        event.save(event)
    end
    println("finished")
    return

end


function demo1()
    option = ""
    inputFile = "cell_test" * option * ".dat"
    inputFile = "new_cell_test.dat"
    println("input:",inputFile)
    abm_data = parseInputFile(inputFile)
    println("================================")
    cellIndex = 0
    outputDir, outputPrefix  = abm_data["output"]
    outPrefix = "demo1"
    println("output:",outputDir, " ", outputPrefix)
    cellCountFile = outputDir * "/cell_count.csv"
    cellCountFile = "cell_count.csv"
    plotFile = outputDir * "/" * outPrefix * "_cell_count.png"
    run__(abm_data,"demo1_",cellIndex)
    println("count:",cellCountFile)
    println("plot:",plotFile)
    plotCellCounts(cellCountFile,plotFile)
    exit()
    return
end

function demo2(growthLimit,inputFile)
    abm_data = parseInputFile(inputFile)
    nCycles = abm_data["nCycles"]
    outputDir, outputPrefix  = abm_data["output"]
    cellIndex = 0
    run__(abm_data,"demo",cellIndex)

    outPrefix = "demo2"
    cellCountFile = "demo_model_count_B.csv"
    targetFile = outputDir * "/" * cellCountFile
    if(isfile(targetFile))
        rm(cellCountFile)
    end
    mv(cellCountFile,targetFile, force=true)
    println(cellCountFile," ",targetFile)
    plotFile = outputDir * "/" * outPrefix * "_cell_count.png"
    cellCountFile = "demo_model_count_A.csv"
    plotCellCounts(cellCountFile,plotFile)
    targetFile = outputDir * "/" * cellCountFile
    if(isfile(targetFile))
        rm(targetFile)
    end
    println("cellCountFile:",cellCountFile)
    println("plotFile:",plotFile)
    plotCellCounts(cellCountFile,plotFile)
    mv(cellCountFile,targetFile, force=true)
    new_dir = "demo2_" * growthLimit * "_" * string(nCycles)
    #createDirectory(newDir)
    mv(outputDir,new_dir)
    println("demo2 finished")
    exit()
    return
end

function renameOutput(outputDir,newOutputDir)

end

function getGrowthLimit(inputFile)
    list = readTextFile(inputFile) 
    index = findText("growthLimit=",list,1)
    index = findText("Events:A",list,1)
    println("index=",index)
    line = strip(list[index])
    index = findfirst("growthLimit=",line)
    text = line[index[1]:end]
    index = findfirst(")",text)
    text = text[1:index[1]-1]
    text = replace(text,"growthLimit=" => "")
    text = strip(text)
    index = findfirst(",",text)
    if(index != nothing)
        text = strip(text[1:index[1]-1])
    end
    return(strip(text))
end

function run_demo1(inputFile)
    Random.seed!(6549)
    println("run_demo1")
    abm_data = parseInputFile(inputFile)
    nCycles = abm_data["nCycles"]
    outputDir, outputPrefix  = abm_data["output"]
    cellIndex = 0
    testIntegrate__(abm_data,"demo",cellIndex)
    println("finished")
    return
end

function run_demo2(inputFile)
    println("run_demo2")
    growthLimit  =  getGrowthLimit(inputFile)
    abm_data = parseInputFile(inputFile)
    nCycles = abm_data["nCycles"]
    outputDir, outputPrefix  = abm_data["output"]
    cellIndex = 0
    run__(abm_data,"demo",cellIndex)

    outPrefix = "demo2"
    cellCountFile = "demo_model_count_B.csv"
    targetFile = outputDir * "/" * cellCountFile
    if(isfile(cellCountFile))
        rm(cellCountFile)
    end
    #println("mv ", cellCountFile, " ", targetFile)
    #mv(cellCountFile,targetFile, force=true)
    plotFile = outputDir * "/" * outPrefix * "_cell_count.png"
    cellCountFile = "demo_model_count_A.csv"
    plotCellCounts(cellCountFile,plotFile)
    println("cellCountFile:",cellCountFile)
    println("plotFile:",plotFile)
    plotCellCounts(cellCountFile,plotFile)
    targetFile = outputDir * "/" * cellCountFile
    if(isfile(targetFile))
        rm(targetFile)
    end
    mv(cellCountFile,targetFile, force=true)
    new_dir = "demo2_" * growthLimit * "_" * string(nCycles)
    mv(outputDir,new_dir)
    println("finished")
    return
end

function parseCommandLine()

    ENV["GKSwstype"] = 100
    #parameterDictionary = CommandParser.readParameters("runCell.dat")
    parameterDictionary = CommandParser.readParameters("options.dat")

    dict = CommandParser.parse(parameterDictionary,ARGS)
    #println(dict)
    inputFile = nothing
    outputDirectory = nothing
    cpus = 1
    demo = nothing
    if (haskey(dict,"-demo"))
       demo = dict["-demo"]
    end
    if (haskey(dict,"-cpus"))
       cpus = dict["-cpus"]
    end
    if (haskey(dict,"-i"))
       inputFile = dict["-i"]
    end
    if (haskey(dict,"-o"))
       outputDirectory = dict["-o"]
    end
    if(demo == "demo1")
        run_demo1("demo1.dat")
    elseif(demo == "demo2")
        run_demo2("demo2.dat")
    else
        println("Error, failed to select a valid demo: {demo1, demo2}")
        exit()
    end
    return
end

function main()
    println("main")
    demo1()
    inputFile = "demo2.dat"
    growthLimit = "0.01"
    demo2(growthLimit,inputFile)
    return
end

#main()

parseCommandLine()

