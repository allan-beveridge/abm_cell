
function do_nothing(model,cell,event)
    return
end

function do_nothing(event)
    return
end

function always_true(model,cell,event)
    return(true)
end

#function dummy__(event)
#end

function modelDummy__(event)
end

function Dummy__(model,cell,event)
    return(true)
end

function dummySave(model,cell,event)
end

mutable struct Event
    name_::String
    data_::Dict{String,Any}
    functions_::Dict{Symbol,Any}
    global_::Bool
end

Base.getproperty(x::Event, property::Symbol) = getfield(x, :functions_)[property]
Base.setproperty!(x::Event, property::Symbol, value) = getfield(x, :functions_)[property] = value
Base.propertynames(x::Event) = keys(getfield(Event, :functions_))

function saveCellEvents(model)
    for cell in allagents(model)
        eventList = getfield(cell,:events_)
        for event in eventList
            event.save(model,cell,event)
        end
    end
end

function updateEvent__(model,cell,event)
    return(true)
end

function update_nutrients(model,cell,event)
    eventData = getfield(event,:data_)
    nutrients = eventData["nutrients"]
    #println("update_nutrients:",model.nutrients)
    for (key,data) in nutrients
        nutrient = data[1]
        value = data[2]
        symbol = Symbol(data[3])
        for cell in allagents(model)
            cellValue = getproperty(cell,symbol)
            value -= cellValue
        end
        nutrients[key] = (nutrient,value,symbol)
    end
    return

end

function countCells(model,cell,event)
    #println("countCells:",model.nSteps)
    cellLineage = getfield(cell,:lineage_)
    #println(cellLineage.status_)
    eventData = getfield(event,:data_)
    agent_id = getfield(cell,:id)
    label = getfield(cell,:label_)
    currentTimeStep, varType = eventData["currentTimeStep"]
    previousTimeStep, varType = eventData["previousTimeStep"]

    if(cellLineage.status_ == "divided")
        timeStepList, varType = eventData["timeStepList" * label]
        cellCountList, varType = eventData["cellCountList" * label]
        cellCount, varType = eventData["cellCount" * label]
        cellCount += 1
        #println(model.nSteps," ==cell:",agent_id, "-", label, " " , cellCount)
        eventData["cellCount" * label] = (cellCount,varType)
        if(model.nSteps == currentTimeStep)
            currentIndex = length(cellCountList)
            cellCountList[currentIndex] = cellCount
        else
            previousTimeStep = currentTimeStep
            currentTimeStep = model.nSteps
            push!(cellCountList,cellCount)
            push!(timeStepList,currentTimeStep)
        end
    end
    return
end

function old_countCells(model,cell,event)
    eventData = getfield(event,:data_)
    agent_id = getfield(cell,:id)
    label = getfield(cell,:label_)
    #println("cell:",agent_id, " ", label)
    indexA = eventData["indexA"]
    indexB = eventData["indexB"]
    previousCellCountA = eventData["previousCellCountA"]
    previousCellCountB = eventData["previousCellCountB"]
    cellCountA = eventData["cellCountA"]
    cellCountB = eventData["cellCountB"]
    currentTimeStep = eventData["currentTimeStep"]
    previousTimeStep = eventData["previousTimeStep"]
    timeStepListA = eventData["timeStepListA"]
    timeStepListB = eventData["timeStepListB"]
    cellCountListA = eventData["cellCountListA"]
    cellCountListB = eventData["cellCountListB"]
    #println("T:",model.nSteps, " ", currentTimeStep)
    if(model.nSteps == currentTimeStep)
        if(label == "A")
            cellCountA += 1
            eventData["cellCountA"] = cellCountA
        elseif(label == "B")
            cellCountB += 1
            eventData["cellCountB"] = cellCountB
        end
    else
        #println("xxx ",previousTimeStep, " ", model.nSteps)
        previousTimeStep = currentTimeStep
        currentTimeStep = model.nSteps
        eventData["currentTimeStep"] = currentTimeStep
        eventData["previousTimeStep"] = previousTimeStep

        if(cellCountA > previousCellCountA)
            #println(indexA, " A count:", cellCountA)
            #println(indexB, " B count:", cellCountB)
            indexA += 1
            eventData["indexA"] = indexA
            push!(cellCountListA,cellCountA)
            push!(timeStepListA,currentTimeStep)
            eventData["previousCellCountA"] = cellCountA
    
            #println(currentTimeStep)
            
            #exit(88)
        end
        if(cellCountB > previousCellCountB)
            indexB += 1
            eventData["indexB"] = indexB
            push!(cellCountListB,cellCountB)
            push!(timeStepListB,currentTimeStep)
            eventData["previousCellCountB"] = cellCountB
        end
    end
    return
    cellCountList = data["cellCountA"]
    timeStepList = data["timeStepA"]
    index = length(cellCountList)
    lastUpdate = timeStepList[index]
    if(model.nSteps != lastUpdate)
        nCells = cellCountList[index] + 1
        push!(cellCountList,nCells)
        push!(timeStepList,model.nSteps)
    else
        exit(88)
    #println("inc count cell:",cell.id)
    end
    #println(cellCountList)
    #println(timeStepList)
    return
end

function saveCellCounts(model,cell,event)
    #println("saveCellCounts")
    eventData = getfield(event,:data_)
    agent_id = getfield(cell,:id)
    label = getfield(cell,:label_)
    timeStepList, varType = eventData["timeStepList" * label]
    cellCountList, varType = eventData["cellCountList" * label]
    df = DataFrame(timeStep=timeStepList,cellCount=cellCountList)
    csvFile = "cell_count_" * label * ".csv"
    CSV.write(csvFile,df)
    return
end

function demoModelCellCount(model,cell,event)
    #println("demoModelCellCount:",model.nSteps)
    countA = 0
    countB = 0
    for cell in allagents(model)
        agent_id = getfield(cell,:id)
        label = getfield(cell,:label_)
        #println("cell:",agent_id, " ", label)
        if(label == "A")
            countA += 1
        elseif(label == "B")
            countB += 1
        end
    end
    eventData = getfield(event,:data_)
    cellCountList, varType = eventData["cellCountA"]
    #println(cellCountList)
    timeStepList, varType = eventData["timeStepA"]
    index = length(cellCountList)
    nCells = cellCountList[index]
    if(countA > nCells)
        push!(cellCountList,countA)
        push!(timeStepList,model.nSteps)
    end
    #println("A:",cellCountList)
    #println("t:",timeStepList)

    cellCountList, varType = eventData["cellCountB"]
    timeStepList, varType = eventData["timeStepB"]
    index = length(cellCountList)
    nCells = cellCountList[index]
    if(countB > nCells)
        push!(cellCountList,countB)
        push!(timeStepList,model.nSteps)
    end
    #println("B:",cellCountList)
    #println("t:",timeStepList)
    return
end

function saveModelCellCount(event)
    #println("saveModelCellCount")
    data = getfield(event,:data_)
    #println(data["timeStepB"])
    #println(data["cellCountB"])
    timeStepList, varType = data["timeStepA"]
    cellCountList, varType = data["cellCountA"]
    df = DataFrame(timeStep=timeStepList,cellCount=cellCountList)
    CSV.write("demo_model_count_A.csv",df)
    timeStepList, varType = data["timeStepB"]
    cellCountList, varType = data["cellCountB"]
    df = DataFrame(timeStep=timeStepList,cellCount=cellCountList)
    CSV.write("demo_model_count_B.csv",df)
    return
end

function saveCellCount(event)
    #println("saveCellCount")
    data = getfield(event,:data_)
    cellCountFile = data["cellCountFile"]
    #println("cellCountFile:",cellCountFile)
    df = DataFrame(timeStep=data["timeStep"],cellCount=data["cellCount"])
    CSV.write(cellCountFile,df)
end

function updateCellCount(model,cell,event)
    nCells = length(model.agents)
    #println("updateCellCount:",nCells, " ", model.nSteps)
    data = getfield(event,:data_)
    cellCountList = data["cellCount"]
    timeStepList = data["timeStep"]
    index = length(cellCountList)
    lastCount = cellCountList[index]
    if(nCells > lastCount)
        push!(cellCountList,nCells)
        push!(timeStepList,model.nSteps)
    end
end

function getModelEvents(nutrients)
    list = Vector{Event}()
    map = Dict{String,Any}()
    #println("n:",nutrients)
# n:Dict{String, Any}("s" => ("s", 1000.0, "s_i", "nutrients:{s=1000.0->s_i}"))
    data = Dict{String,Any}()
    data["nutrients"] = nutrients
    functions = Dict{Symbol,Any}()
    functions[:test] = updateEvent__
    functions[:execute] = update_nutrients
    functions[:save] = modelDummy__
    event = Event("updateNutrients",data,functions,false)
    push!(list,event)

    data = Dict{String,Any}()
    data["cellCount"] = [0]
    data["timeStep"] = [0]
    functions = Dict{Symbol,Function}()
    functions[:test] = updateEvent__
    functions[:execute] = updateCellCount
    functions[:save] = saveCellCount
    event = Event("updateCellCount",data,functions,false)
    push!(list,event)
    return(list)
end

function getCellEvents()
    println("getCellEvents")
    list = Vector{Event}()
    data = Dict{String,Any}()
    data["maxSteps"] = (5, Int64)
    data["maxSteps"] = (50, Int64)
    data["fraction"] = (0.5, Int64)
    data["nutrient"] = (0.00005, Float64)
    functions = Dict{Symbol,Any}()
    testFunction = convertToFunction("timeToDivide")
    functions[:test] = timeToDivide
    functions[:execute] = divide_cell
    functions[:save] = modelDummy__
    event = Event("divideByTime",data,functions,false)
    push!(list,event)
    return(list)
end

function timeToDivide(model,cell,event)
    #println("timeToDivide:",event.data_, " maxSteps:",event.data_["maxSteps"])
    eventData = getfield(event,:data_)
    #println("d:",eventData)
    nutrient, varType = eventData["nutrient"]
    nCells = length(model.agents)
    #max = 0.05 - nCells * 0.00005
    max = 0.05 - nCells * nutrient
    nSteps = getfield(cell,:nSteps_)
    probability = nSteps * rand(Uniform(0.0,max))
    if(probability > 1.0)
        #println("timeToDivide:",model.nSteps)
        return(true)
    end
    return(false)
end

function divide_cell(model,cell,event)
    agent_id = getfield(cell,:id)
    #println("DIVIDE CELL:",agent_id, " ", typeof(agent_id))
    nextAgent = model.lastCell + 1
    p_id = getproperty(cell,:id)
    cell_id = getproperty(cell,:id_)
    cellIndex = cell_id - agent_id
    model.lastCell += 2

    name = getfield(event,:name_)
    data = getfield(event,:data_)
    fraction, varType = data["fraction"]
    u1 = divideResources(cell,fraction)
    u2 = divideResources(cell,1.0-fraction)

    cellLineage = getfield(cell,:lineage_)
    cellLineage.status_ = "divided"
    cell1 = createNewCell(model,cell,nextAgent,u1,(0.0,0.0),cellIndex)
    cell2 = createNewCell(model,cell,nextAgent+1,u2,(0.0,0.0),cellIndex)
    kill_agent!(cell,model)
    return()
end

    function divideResources(cell,fraction)
        nSymbols = length(cell.u0_)
        u = zeros(nSymbols)
        for i in 1:nSymbols
            u[i] = fraction * getproperty(cell,cell.u0_[i])
        end
        return(u)
    end

    function convertToPairedList_(cell,u0List)
        integrator = getfield(cell,:integrator_)
        u0_vec = Vector{Pair{Symbol, Float64}}()
        index = 1
        for field in cell.u0_
            push!(u0_vec,Pair(field,u0List[index]))
            index += 1
        end
        return(u0_vec)
    end

    function createNewCell(model,cell,agent_id,u0,pos,cellIndex)
        #println("createNewCell:",agent_id, " ", cellIndex)
        cell_id = agent_id + cellIndex

        #u0 = convertPairsToFloats(cell.u0_)
        _rn = model.rn
        _p = model.p
        _u0 = model.u0

        _rn = cell.rn_
        _p = cell.p_
        _u0 = convertToPairedList_(cell,u0)
#println("_u0:",_u0)

        prob = ODEProblem(_rn,u0,model.tspan,_p)
        integrator = OrdinaryDiffEq.init(prob, OrdinaryDiffEq.Rodas4(); advance_to_tstop = true)
        d = getCellDictionary(agent_id,cell_id,pos,_u0,integrator)

        cellLineage = getfield(cell,:lineage_)
        ancestors = deepcopy(cellLineage.ancestors_)
        push!(ancestors,string(cell.id))
        created = cellLineage.timeSteps_
        lineage = Lineage(string(agent_id), created,cellLineage.deltaT_,"alive",1,ancestors)

        events = copyCellEvents(cell)
        label = getfield(cell,:label_)
        newCell = add_agent!(pos,model,d,lineage,integrator,events,label,0)
        newCell.p_ = deepcopy(cell.p_)
        newCell.rn_ = deepcopy(cell.rn_)
        return(newCell)
    end

function copyCellEvents(cell)
    #println("copyEvents")
    cellEvents = getfield(cell,:events_)
    #println(cellEvents)
    eventList = Vector{Any}()
    label = getfield(cell,:label_)
    #println("Lbel:",label)
    for event in cellEvents
        _global = getfield(event,:global_)
        #println("g:",_global)
        if(_global)
            copy = event
        else
            copy = deepcopy(event)
        end
        push!(eventList,copy)
    end
    return(eventList)
end

function updateCellCount(model,cell,event)
    nCells = length(model.agents)
    data = getfield(event,:data_)
    cellCountList = data["cellCount"]
    timeStepList = data["timeStep"]
    index = length(cellCountList)
    lastCount = cellCountList[index]
    if(nCells > lastCount)
        push!(cellCountList,nCells)
        push!(timeStepList,model.nSteps)
    end
end

module Events
    include("parser.jl")

struct Equation
    formula_::String
    variables_::Vector{String}
    symbolChar::String
end

function parseEvents(eventMap::Dict{String,Main.Event},textList)
    nLines = length(textList)
    index = 0
    eventIndex = 0
    eventPrefix = nothing
    while(index < nLines)
        index += 1
        line = textList[index]
        text = strip(line)
        if(startswith(text,"ModelEvent:"))
            eventIndex = index
            eventPrefix = "ModelEvent:"
            break
        elseif(startswith(text,"CellEvent:"))
            eventPrefix = "CellEvent:"
            eventIndex = index
            break
        end
    end
    if(eventIndex == 0)
        return("no events found")
    end
    while index < nLines
        index += 1
        line = textList[index]
        text = strip(line)
        if(startswith(text,"ModelEvent:") )
            event = parseEvent_(textList[eventIndex:index-1])
            eventIndex = index
            name = getfield(event,:name_)
            eventMap[getfield(event,:name_)] = event
        elseif(startswith(text,"CellEvent:") )
            event = parseEvent_(textList[eventIndex:index-1])
            eventIndex = index
            name = getfield(event,:name_)
            eventMap[getfield(event,:name_)] = event
        elseif(startswith(text,"end"))
            event = parseEvent_(textList[eventIndex:index-1])
            eventMap[getfield(event,:name_)] = event
            break
        end
    end
    return(nothing)
end

function parseEventFunction__(prefix,text)
    functionName = replace(text,prefix => "")
    functionName = strip(functionName)
    try
        fnc = getfield(Main, Symbol(functionName))
        #println("fnc:",fnc)
        return(fnc)
    catch e # UndefVarError
        println("parseEventFunction__:",text)
        error = "Invalid function:" * functionName * "\n" * sprint(showerror,e)
        abort__(error * "\n")
    end
end

function parseEvent_(textList)
    #println("parseEvent_")
    text = strip(textList[1])
    cols = split(strip(text),":")
    if(length(cols) != 2)
       println("Error:",textList[1])
       abort__(textList[1],"Error, invalid 'Event'")
    end
    label = strip(cols[2])
    testFunction = getfield(Main, Symbol("always_true"))
    executeFunction = nothing
    saveFunction = getfield(Main, Symbol("do_nothing"))
    global_ = false
    map = Dict{String,Any}()
    nLines = length(textList)
    for i in 2:nLines
        line = textList[i]
        text = strip(line)
        if(startswith(text,"data:"))
            cols = split(text,":")
            if(length(cols) != 2)
                abort__(line,"Error, invalid 'Event' data")
            end
            var_cols = split(cols[2],";")
            for col in var_cols
                varData = parseVariable(col)
                if(varData == nothing)
                    abort__(line,"Error, invalid 'Event' variable")
                end
                map[varData[1]] = (varData[2],varData[3])
            end
        elseif(startswith(text,"equation:"))
            eqnLabel,eqn = parseEquation_(text)
            map[eqnLabel] = eqn
        elseif(startswith(text,"test:"))
            testFunction = parseEventFunction__("test:",text)
        elseif(startswith(text,"execute:"))
            executeFunction = parseEventFunction__("execute:",text)
        elseif(startswith(text,"save:"))
            saveFunction = parseEventFunction__("save:",text)
        elseif(text == "copy")
            global_ = false
        elseif(text == "global")
            global_ = true
        else
            abort__(line,"Error, invalid 'Event' input")
        end
    end
    functions = Dict{Symbol,Any}()
    functions[:test] = testFunction
    functions[:execute] = executeFunction
    functions[:save] = saveFunction
    event = Main.Event(label,map,functions,global_)
    return(event)

end 

function parseEquation_(text)
    text = replace(text, "equation:" => "")
    cols = split(strip(text),"=")
    if(length(cols) != 2)
        abort__("\nError parsing equation:" * text * "\n")
    end
    label = strip(cols[1])
    equation = strip(cols[2])
    eqn = Events.generateEquation(equation)
    return(label,eqn)
end

function splitString(s,substr)
    index = findfirst(substr,s)
    #posn = index[1] + length(substr) - 1
    #posn = index[end]
    s1 = s[1:index[1]-1]
    s2 = s[index[1]:end]
    s2 = replace(s2,substr => "",count=1)
    #s2 = s[posn+1:end]
    return(s1,s2)
end

function generateEquation(equation,symbolChar = "\$")
    println("generateEqn:",equation)
    char = " "
    eqn = replace(equation,"(" => " ")
    eqn = replace(eqn,")" => char)
    eqn = replace(eqn,"+" => char)
    eqn = replace(eqn,"-" => char)
    eqn = replace(eqn,"*" => char)
    eqn = replace(eqn,"/" => char)
    cols = split(strip(eqn))
    start = 1
    eqn = ""
    s2 = equation * " "
    posn = 1
    for key in cols
        s1, s2 = splitString(s2,key)
        eqn *= s1 * symbolChar * key
    end
    return(Events.Equation(eqn,cols,symbolChar))
end

function evaluateEquation(equation,parameterMap,symbolChar = "\$")
    #println("p:",parameterMap)
    #println("e:",equation)
    eqn = equation.formula_
    for var in equation.variables_
        key = symbolChar * var
        value = parameterMap[key]
        #println(var, " ", value, " ", typeof(value))
        eqn = replace(eqn,key => value)
    end
    #println(eqn)
    ex = Meta.parse(eqn)
    res = eval(ex)
    return(res)
end

function evaluate__(equation,variableMap,symbolChar="\$")
    println("evaluate__:",equation)
    eqn = equation.formula_
    for key in equation.variables_
        if (!haskey(variableMap,key))
            abort__("missing variable:" * key * "\n")
        end
        value = variableMap[key]
        #println(key,":",value)
        eqn = replace(eqn,symbolChar * key => value)
    end
    println(equation.formula_)
    println("evaluate eqn:",eqn)
    try
        ex = Meta.parse(eqn,raise=true)
        result = eval(ex)
        return(result)
    catch e # UndefVarError
        println("Error, Events::evaluate__:",equation.formula_)
        error = "Invalid function:" * functionName * "\n" * sprint(showerror,e)
        abort__(error * "\n")
    end
end

function getCellParameters(cell,symbolChar="\$")
    map = Dict{}()
    for s in cell.u0_
        map[symbolChar*string(s)] = getproperty(cell,s)
    end
    for (key,value) in cell.p_
        map[symbolChar * string(key)] = value
    end
    return(map)
end


end # module

function parseEquation_(text)

    text = replace(text, "equation:" => "")
    println(text)
    cols = split(strip(text),"=")
    println(cols)

end

function test()
    eqn = "λ = ((c_q + c_m + c_t + c_r) * (gmax*a/(Kgamma + a)))/M"
    eqn = "((c_q + c_m + c_t + c_r) * (gmax*a/(Kgamma + a)))/M"
    #eqn = "(c_q + c_m + c_t + c_r) * (gmax*a/(Kgamma + a)))/M"
    map = Dict("c_q" => 1, "c_m" => 1, "c_t" => 1, "c_r" => 1, "gmax" => 1, "kGamma" => 1, "a" => 1, "M" => 1)
    map["Kgamma"] = 1.0
    equation = Events.generateEquation(eqn)
    println(equation)
    value = Events.evaluate__(equation,map)
    println("v:",value)
end

function main()

    #test()
    #exit()

    text = "equation:growthRateEqn = ((c_q + c_m + c_t + c_r) * (gmax*a/(Kgamma + a)))/M"
    label,eqn = Events.parseEquation_(text)
    println(label)
    println(eqn)
end


#main()
