# https://stackoverflow.com/questions/55040847/julia-get-datatype-from-string

cellFunctionFile = "UserCellFunctions.jl"
if(isfile(cellFunctionFile))
    include(cellFunctionFile)
end

macro datatype(str); :($(Symbol(str))); end

mutable struct CellEventData
    name_::String
    data_::Dict{String,Any}
    testFunction::Function
    executeFunction::Function
end

function parseCellEvent(line)
    line = replace(line,"cell_event:" => "")
    line = strip(line)
    line = line[2:length(line)-1]
    columns = split(line,",")
    eventName = strip(columns[1])
    varMap = Dict{String,Any}()
    testFunctionName = parseFunctionData(varMap,columns[2])
    symbol = Symbol(testFunctionName)
    testFunction = getfield(Main, symbol)
    #println("F:",testFunctionName,"::",varMap)
    eventFunctionName = parseFunctionData(varMap,columns[3])
    symbol = Symbol(eventFunctionName)
    eventFunction = getfield(Main, symbol)
    #println("F:",eventFunctionName,"::",varMap)
    event = CellEventData(eventName,varMap,testFunction,eventFunction)
    return(event)
end

function parseFunctionData(varMap,line)
    line = strip(line)
    columns = split(line,"[")
    if(length(columns) != 2)
        println("parseFunctionData error-1")
        exit(999)
    end
    functionName = strip(columns[1])
    if(last(columns[2]) != ']')
        println("parseFunctionData error-2")
        exit(999)
    end
    functionData = columns[2][1:length(columns[2])-1]
    columns = split(functionData,";")
    for col in columns
        parseFunctionVariable(varMap,col)
    end
    return(functionName)
end

function parseFunctionVariable(varMap,data)
    data = strip(data)
    #println("v:",data)
    columns = split(data,"=")
    if(length(columns) != 2)
        println("parseFunctionVariable: error-1")
        exit(999)
    end
    varName = columns[1]
    cols = split(columns[2],"(")
    if(length(cols) != 2)
        println("parseFunctionVariable: error-2")
        exit(999)
    end
    cols[1] = strip(cols[1])
    cols[2] = strip(cols[2])
    if(last(cols[2]) != ')')
        println("parseFunctionVariable: error-3")
        exit(999)
    end
    varType = cols[2][1:length(cols[2])-1]
    type = getfield(Base, Symbol(varType))
    #type = @datatype(varType)
    value = parse(type,cols[1])
    varMap[varName] = (value,type)
end

function DivideByTime(model,cell,event)
    #println("DivideByTime")
    nSteps = getfield(cell,:nSteps_)
    maxSteps, varType = event.data_["maxSteps"]
    if(nSteps > maxSteps)
        return(true)
    end
    return(false)
end

function DivideCell(model,cell,event)
    cell_id = getfield(cell,:id)
    #println("DivideCell:",cell_id, " ", typeof(cell_id))
    nextCell = model.lastCell + 1
    model.lastCell += 2

    fraction, type = event.data_["fraction"]
    u1 = divideResources(cell,fraction)
    u2 = divideResources(cell,1.0-fraction)
    #println("u1:",u1)
    #println("u2:",u2)

    cellLineage = getfield(cell,:lineage_)
    cellLineage.status_ = "divided"
    cell1 = createNewCell(model,cell,nextCell,u1,(0.0,0.0))
    cell2 = createNewCell(model,cell,nextCell+1,u2,(0.0,0.0))
    kill_agent!(cell,model)
    return()
end


#test = "cell_event:{Divide Cell,DivideByTime[maxSteps=5(Int64);x=10(Float64)],DivideCell[fraction=0.5(Float64)]}"
#parseCellEvent(test)
