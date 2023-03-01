using Agents
using DataFrames

#include("parser.jl")
include("parseCell.jl")
include("Cell.jl")

struct abmData
    cellMap_
    cellEventsMap_
    EventsMap_ # stores all defined events
    modelEventsMap_
    modelEventsList_
    cellList_
end

function parseEvents(eventMap::Dict{String,Main.Event},list)
    #println("==parseEvents")
    nLines = length(list)
    index = 0
    while index < nLines
        index += 1
        line = strip(list[index])
        if(startswith(line,"Events"))
            index2 =  findText("end",list,index)
            if(index2 != nothing)
                error = Events.parseEvents(eventMap,list[index:index2])
                return(error)
            else
                error = "parseEvents() error, invalid input:" * line * "\n"
                return(error)
            end
            break
        end
    end
    return(nothing)
end

function readCells(cellFile)
    println("reading:",cellFile)
    if(!isfile(cellFile))
        abort__("Error, invalid file:" * cellFile * "\n")
    end
    list = readTextFile(cellFile)
    list = removeBlankLines(list)
    list = removeCommentLines(list,"#")
    return(list)
end

function parseCells(textList)
    cellList = Vector{Any}()
    for line in textList                # Cell A 1 (0.0,0.0,0.0)
        cols = split(strip(line))
        cellType = strip(cols[2])
        cell_id = Base.parse(Int64,strip(cols[3]))
        coords = strip(cols[4])
        coords = replace(coords,"(" => "")
        coords = replace(coords,")" => "")
        cols = split(coords,",")
        x = Base.parse(Float64,strip(cols[1]))
        y = Base.parse(Float64,strip(cols[2]))
        #println(cellType," ", cell_id, " ", x, " ", y)
        data = (cellType,cell_id, (x,y))
        push!(cellList,data)
    end
    return(cellList)
end

function readModelEventsFile(eventsFile)
    println("reading:",eventsFile)
    if(!isfile(eventsFile))
        abort__("Error, invalid file:" * eventsFile * "\n")
    end
    list = readTextFile(eventsFile)
    list = removeBlankLines(list)
    list = removeCommentLines(list,"#")
    eventMap = Dict{String,Main.Event}()
    parseEvents(eventMap,list)
    return(eventMap)
end

function readEventsFile(eventsMap::Dict{String,Main.Event},eventsFile)
    println("reading:",eventsFile)
    if(!isfile(eventsFile))
        abort__("Error, invalid file:" * eventsFile * "\n")
    end
    list = readTextFile(eventsFile)
    list = removeBlankLines(list)
    list = removeCommentLines(list,"#")
    error = parseEvents(eventsMap,list)
    if(error != nothing)
        message = "Error reading:" * eventsFile * "\n" * error * "\n"
        abort__(message)
    end
    return(nothing)
end

function readCellEventsFile(eventsFile)
    if(!isfile(eventsFile))
        abort__("Error, invalid file:" * eventsFile * "\n")
    end
    list = readTextFile(eventsFile)
    list = removeBlankLines(list)
    list = removeCommentLines(list,"#")
    eventMap = parseEvents(eventsMap::Dict{String,Main.Event},list)
    println(eventMap)
    return(eventMap)
end

function parseCellEvents(inputLine,columnList)
    #println("parseCellEvents:",columnList)
    eventList = Vector{}()
    for i in 1:length(columnList)
        eventData = parseCellEvent__(inputLine,columnList[i])
        push!(eventList,eventData)
    end
    return(eventList)
end

function parseCellEvent__(inputLine,text)
    cols = split(text,"{")
    if(length(cols) != 2)
        abort__(inputLine,"Error, invalid input:" * text * "\n")
    end
    if( !endswith(cols[2],"}") )
        abort__(inputLine,"Error, invalid input:" * text * "\n")
    end
    eventLabel = strip(cols[1])
    cols = split(strip(cols[2]),"}")
    if(length(cols) != 2)
        abort__(inputLine,"Error, invalid input:" * text * "\n")
    end
    eventName,global_,variableList = parseEventData__(inputLine,strip(cols[1]))
    return(eventLabel,eventName,global_,variableList)
end

function parseEventData__(inputLine,data)
    text = strip(data)
    #println("parseEventData__::",text)
    variableList = Vector{Pair{String,String}}()
    global_ = false
    if(endswith(text,")") )
        text = text[1:length(text)-1]
        cols = split(text,"(")
        if(length(cols) != 2)
            abort__(inputLine,"Error, invalid input:" * text * "\n")
        end
        eventName = strip(cols[1])
        columns = split(strip(cols[2]),",")
        for col in columns
            cols = split(col,"=")
            if(length(cols) > 2)
                abort__(inputLine,"Error, invalid input:" * col * "\n")
            elseif(length(cols) == 2)
                pair = Pair(strip(cols[1]),strip(cols[2]))
                push!(variableList,pair)
                #variableMap[strip(cols[1])] = strip(cols[2])
            elseif(strip(cols[1]) == "global")
                global_ = true
            end
        end
    else
        eventName = text
    end
    return(eventName,global_,variableList)
end

function old_parseCellEvents(cols)
    eventList = Vector{String}()
    for i in 1:length(cols)
        eventName = strip(cols[i])
        push!(eventList,strip(cols[i]))
    end
    return(eventList)
end

function parseDictionary(dictionaryFile)
    println("reading:",dictionaryFile)
    if(!isfile(dictionaryFile))
        abort__("Error, invalid file:" * dictionaryFile * "\n")
    end
    list = readTextFile(dictionaryFile)
    list = removeBlankLines(list)
    list = removeCommentLines(list,"#")
    map = Dict{String,Any}()
    for line in list
        line = strip(line)
        if(!endswith(line,")"))
            abort__("Invalid input:" * line * "\n")
        end
        name, dataType = getVariableType(line)
        if(dataType == nothing)
            abort__("Invalid input:" * line * "\n")
        end
        map[name] = dataType
    end
    return(map)
end

function getVariableType(text)
    text = strip(text)
    if( endswith(text,")") )
        text = text[1:length(text)-1]
        cols = split(text,"(")
        if(length(cols) != 2)
            write(Base.stderr,"getVariableType::Error, invalid var data:")
            return(nothing,nothing)
        end
        dataType =  getDataType(cols[2])
        return(strip(cols[1]),dataType)
    else
        println("getVariableType::Error, invalid var data:",text)
        exit(999)
    end
    return(nothing)
end

function getDataType(varType)
    #println("getDataType:",varType)
    code = "dataType = " * strip(varType)
    try
        ex = Meta.parse(code)
        dataType = eval(ex)
        return(dataType)
    catch e # UndefVarError
        error = sprint(showerror,e)
        write(Base.stderr,error)
        return(nothing)
    end
end

function parseInputFile(inputFile)
    if(!isfile(inputFile))
        abort__("Error, invalid file:" * inputFile)
    end
    nCycles = 200
    #cellEventsMap = Dict{String,Vector{String}}()
    cellEventsMap = Dict{String,Any}()
    #cellEventsMap = Dict{String,(String,Bool,Vector{Pair{String, String}})}()
    cellEventDataList = []
    modelEventDataList = []
    output = nothing
    modelEventsMap = nothing
    modelEvents = nothing
    cellList = nothing
    steadyStateList = nothing
    #modelEventsList = nothing
    EventsMap = Dict{String,Main.Event}()
    modelEventsMap = Dict{String,Main.Event}()
    cellMap = Dict{String,Any}()
    u0_map = Dict{String,Vector{Pair{Symbol,Float64}}}()
    list = readTextFile(inputFile)
    list = removeBlankLines(list)
    list = removeCommentLines(list,"#")
    dataMap = Dict{String,Any}()
    for line in list
        text = strip(line)
println(text)
        text = removeComment(strip(line),"#")
        cols = split(strip(text),":")
        if(length(cols) != 2)
            abort__(text,"Error, invalid input\n")
        end
        command = strip(cols[1])
        if(command == "Cell")
            cols = split( strip(cols[2]) )
            if(length(cols) != 2)
                abort__(text,"Error, invalid input\n")
            end
            cellType = strip(cols[1])
            cellFile = strip(cols[2])
            if( endswith(cellFile,".dat") )
                modelParams = parseCell(cellFile)
                cellMap[cellType] = modelParams
            elseif(endswith(cellFile,".dict"))
                cols = split(strip(line))
                label = strip(cols[1])
                label = replace(label,"Cell:" => "")
                dictionaryFile = strip(cols[2])
                map = parseDictionary(dictionaryFile)
                dataMap[label] = map
            end
        elseif(command == "ModelEvents")
            eventsFile = checkFilename__(text,cols[2])
            #modelEventsMap = readModelEventsFile(eventsFile)
            readEventsFile(modelEventsMap,eventsFile)
        elseif(command == "CellEvents")
            eventsFile = checkFilename__(text,cols[2])
            readEventsFile(EventsMap,eventsFile)
        elseif(command == "Events")
            columns = split(cols[2])
            type = strip(columns[1])
            eventList = parseCellEvents(text,columns[2:length(columns)])
            println(eventList)
            if(type == "model")
                modelEventDataList = vcat(modelEventDataList,eventList)
            else
                cellEventsMap[type] = eventList
            end
        elseif(command == "Initialise")
            cellList = readCells(strip(cols[2]))
        elseif(command == "Output")
            output = parseOutput(line)
            if(output[1] == nothing)
                abort__(line,"Error, invalid input")
            end
        #elseif(startswith(line,"Dictionary:"))
        #    cols = split(strip(line))
        #    label = strip(cols[1])
        #    dictionaryFile = strip(cols[2])
        #    map = parseDictionary(dictionaryFile)
        #    println(label,":",map)
        #    exit()
        elseif(command == "u0")
            columns = split(strip(cols[2]))
            cellType = columns[1]
            u0 = parse_u0__(line,columns[2:end])
            u0_map[cellType] = u0
            println(cellType,":",u0)
        elseif(command == "nCycles")
            nCycles = parseValue(Int64,strip(cols[2]))
            if(nCycles == nothing)
                abort__(line,"Error, expected integer:" * cols[2])
            end
            println("nCycles:",nCycles)
        elseif(command == "ENV")
            println(cols)
            key = cols[1]
            value = cols[2]
            println(key, " ", value)
        else
            abort__(line,"parseInputFile() Error, invalid key word:" * command * "\n")
        end
    end
    if(output == nothing)
        abort__("Error, output has not been specified")
    end
    checkEvents__(modelEventDataList,modelEventsMap)
    println(cellEventsMap)
    println(keys(cellEventsMap))
    for (key,eventDataList) in cellEventsMap
        println(key,":",eventDataList)
        checkEvents__(eventDataList,EventsMap)
    end

    dataMap["u0"] = u0_map
    dataMap["nCycles"] = nCycles
    dataMap["output"] = output
    dataMap["cellMap"] = cellMap
    dataMap["cellEventsMap"] = cellEventsMap
    dataMap["EventsMap"] = EventsMap
    dataMap["modelEventsMap"] = modelEventsMap
    #dataMap["modelEventsList"] = modelEventsList
    dataMap["modelEventsList"] = modelEventDataList
    dataMap["cellList"] = cellList
    return(dataMap)
end

function checkEvents__(eventDataList,eventMap)
    for eventData in eventDataList
        println(eventData)
        eventName = eventData[2]
        if(!haskey(eventMap,eventName))
            println("em:",eventMap)
            abort__("Error, invalid (duplicate) event:" * eventName * "\n")
        end
    end
    return(nothing)
end

function checkFilename__(inputLine,filename)
    filename = strip(filename)
    cols = split(filename)
    if(length(cols) != 1)
       abort__(inputLine,"Error, invalid filename:" * filename * "\n")
    end
    if(!isfile(filename))
       abort__(inputLine,"Error, invalid file:" * filename * "\n")
    end
    return(filename)
end

function parse_u0__(line,columns)
    println("parse_u0:",columns)
    nCols = length(columns)
    scalingFactor = 1.0
    csvFile = columns[1]
    if(nCols < 3)
        checkFilename__(line,csvFile)
        if(nCols == 2)
            try
                scalingFactor = Base.parse(Float64,columns[2])
            catch e # UndefVarError
                error = "Invalid float:" * columns[2] * "\n" * sprint(showerror,e)
                abort__(error * "\n")
            end
        end
    else
       abort__(line,"Error, invalid input:" * columns[2] * "\n")
    end
    pairList = readSteadyState(csvFile,scalingFactor)
    return(pairList)
end

function parseValue(dataType,text)
    try
        value = Base.parse(dataType,text)
        return(value)
    catch e # ArgumentError
        write(Base.stderr,sprint(showerror,e) * "\n")
        return(nothing)
    end
end

function parseOutput(text)
    text = replace(text,"Output:" => "")
    cols = split(text)
    if(length(cols) != 2)
        return(nothing,nothing)
    end
    output = (strip(cols[1]),strip(cols[2]))
    return(output)
end

function convert_to_abm_data(dataMap)
    cellMap =  dataMap["cellMap"]
    cellEventsMap = dataMap["cellEventsMap"]
    EventsMap = dataMap["EventsMap"]
    modelEventsMap = dataMap["modelEventsMap"]
    modelEventsList = dataMap["modelEventsList"]
    cellList = dataMap["cellList"]
    abm_data = abmData(cellMap,cellEventsMap,EventsMap,modelEventsMap,modelEventsList,cellList)
    return(abm_data)
end

function getModelEvents(abm_data)
    EventsMap_ = abm_data["EventsMap"]
    cellEventsMap = abm_data["cellEventsMap"]
    modelEventsMap = abm_data["modelEventsMap"]
    modelEventList = abm_data["modelEventsList"]
    #modelEventList = cellEventsMap["model"]
    eventList = Vector{Any}()
    eventsMap = OrderedDict{String,Event}()
    currentEvent = nothing
    for data in modelEventList
        #println("data:",data)
        eventLabel = data[1]
        eventName = data[2]
        global_ = data[3]
        variableList = data[4]
        if(!haskey(modelEventsMap,eventName))
            abort__("Error, invalid event:" * eventName)
        end
        if(haskey(eventsMap,eventName))
            currentEvent = eventsMap[eventName]
        else
            currentEvent = modelEventsMap[eventName]
        end
        if(global_)
            event = currentEvent
        else
            event = deepcopy(currentEvent)
        end
        resetEvent(event,variableList)
        eventsMap[eventLabel] = event
        #println("event:",event)
    end
    return(eventsMap)
end

function old_getModelEvents(abm_data)
    EventsMap_ = abm_data["EventsMap"]
    cellEventsMap = abm_data["cellEventsMap"]
    modelEventsMap = abm_data["modelEventsMap"]
    modelEventList = cellEventsMap["model"]
    eventList = Vector{Any}()
    for label in modelEventList
        println("Label:",label)
        event = modelEventsMap[label]
        push!(eventList,event)
    end
    return(eventList)
end

function getCellEvents(abm_data,cellType)
    eventsMap = abm_data["EventsMap"]
    eventsMapCopy = deepcopy(eventsMap)
    if(eventsMap == nothing)
        write(Base.stderr,"Error, getCellevents()\n")
        write(Base.stderr,"no events defined\n\n")
        exit(999)
    end
    #println("getCellEvents:",cellType)
    cellEventsMap = abm_data["cellEventsMap"]
    cellEventList = cellEventsMap[cellType]

    eventsMap = generateEvents__(cellEventList,eventsMapCopy)
    eventList = collect(values(eventsMap))
    return(eventList)

end

function generateEvents__(eventDataList, originalEventsMap)
    #println("generateEvents__:",eventDataList)
    eventList = Vector{Any}()
    eventsMap = OrderedDict{String,Event}()
    
    currentEvent = nothing
    for data in eventDataList
        println(data)
        eventLabel = data[1]
        eventName = data[2]
        global_ = data[3]
        variableList = data[4]
        #println(eventName, " ", global_, " ", variableList)
        #println("G:",global_)
        if(!haskey(originalEventsMap,eventName))
            abort__("Error, invalid event:" * eventName)
        end
        if(haskey(eventsMap,eventName))
            currentEvent = eventsMap[eventName]
        else
            currentEvent = originalEventsMap[eventName]
        end
        if(global_)
            event = currentEvent
        else
            event = deepcopy(currentEvent)
        end
        resetEvent(event,variableList)
        eventsMap[eventLabel] = event
        #println("event:",event)
    end
    return(eventsMap)
end

function resetEvent(event,variableList)
    eventData = getfield(event,:data_)
    #println(event)
    #println("data:",eventData)
    for pair in variableList
        if (haskey(eventData,pair.first))
            currentValue, varType = eventData[pair.first]
            dataType = getDataType(varType)
            println(pair.first, ":", pair.second, " ", dataType)
            value = Base.parse(dataType,pair.second)
            eventData[pair.first] = (value,dataType)
        else
            println("ERROR")
            abort__("Error, invalid variable:" * pair.first)
        end
    end
    return(nothing)
end

function initialiseModel(abm_data,cellIndex=0)

    println("initialiseModel")
exit()
    println(abm_data)
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
    u0 = convertPairsToFloats(data.u0_)
    println(u0)
    println(typeof(u0))

    fieldList = getSymbolList(data.u0_)
    println(fieldList)
    println(typeof(fieldList))
    uList = nothing
    cellLabel = "label"
    println(keys(abm_data))
    cellList = parseCells(abm_data["cellList"])
    println(cellList)
    for i in 1:model.nAgents
        println(i, " ", cellList[i])
    end
    for i in 1:model.nAgents
        println(cellList[i])
        cellType, cell_id, position = cellList[i]
        events = getCellEvents(abm_data,cellType)
        println(cellType, " ", cell_id, " ", position)
        println("events:",events)

        prob = ODEProblem(model.rn, u0, data.tspan_, data.p_)
        integrator = OrdinaryDiffEq.init(prob, OrdinaryDiffEq.Rodas4(); advance_to_tstop = true)
        lineage = Lineage(string(i),1,data.deltat_,"alive",1,Vector{String}())
        cell_id = cellIndex + i
        d = getCellDictionary(i,cell_id,position,data.u0_,integrator)
        d[Symbol("label__")] = 2
        cell = add_agent!(position,model,d,lineage,integrator,events,cellLabel,0)
    end
    modelEventsMap =  getModelEvents(abm_data)
    modelEventsList = Vector{Event}()
    for (key,event) in modelEventsMap
        push!(modelEventsList,event)
    end
    model.events =  modelEventsList
    return(model)
end

function run(abm_data,cellIndex=0)
    Random.seed!(6549)
    outputDir = "./test_dir"
    model = initialiseModel(abm_data,cellIndex)

    cellMap = abm_data["cellMap"]
    data_a = cellMap["A"]
    modelParams = data_a
    maxTime = modelParams.tspan_[2]
    nCycles = Int64(maxTime/modelParams.deltat_)
    nCycles = 10
    nCycles = 500
    #nCycles = 200
    aList = [:r,:e_t,:e_m,:q,:m_r,:m_t,:m_m,:m_q,:c_r,:c_t,:c_m,:c_q,:a,:s_i]
    aList = [:id_,:r,:e_t,:e_m,:q,:m_r,:m_t,:m_m,:m_q,:c_r,:c_t,:c_m,:c_q,:a,:s_i]
    #aList = [:label__,:r,:e_t,:e_m,:q,:m_r,:m_t,:m_m,:m_q,:c_r,:c_t,:c_m,:c_q,:a,:s_i]
    agentResults, modelResults = run!(model, update_cell!, update_model!, nCycles; mdata = [], adata = aList)
    CSV.write("zzz_abm_cell.csv",agentResults)

    roundColumns(agentResults,aList)
    CSV.write("general_abm_cell.csv",agentResults)
    plotPrefix = "agent_model_dt=" * string(modelParams.deltat_) * "_maxt=" * string(maxTime)
    println("prefix:",plotPrefix)
    agentPlot(agentResults,model.nAgents,outputDir,plotPrefix)

    eventList = model.events
    for event in eventList
        println("event:",event)
        event.save(event)
    end
    println("finished")
    return

end

function readSteadyState(csvFile, scalingFactor=1.00)
    df = DataFrame(CSV.File(csvFile))
    labelList = df[:,:label]
    valueList = df[:,:value]
    nValues = length(valueList)
    pairList = Vector{Pair{Symbol,Float64}}()
    for i in 1:nValues
        pair = Pair(Symbol(labelList[i]),scalingFactor*valueList[i])
        push!(pairList,pair)
    end
    return(pairList)
end

function test()
    inputFile = "cell_test.dat"
    dataMap = parseInputFile(inputFile)
exit()
    abm_data = convert_to_abm_data(dataMap)
    cellIndex = 0
    run(abm_data,cellIndex)
    return
end

function test_parse()
    inputFile = "new_cell_test.dat"
    dataMap = parseInputFile(inputFile)
    return
end

function main()
    #test()
    #test_parse()
    csvFile = "steadystate.csv"
    pairList = readSteadyState(csvFile)
    println(pairList)
end

#main()
