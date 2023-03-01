using Catalyst

include("parser.jl")
include("Events.jl")
include("CellFunctions.jl")

struct model_params
    nCycles_::Int64
    nCells_::Int64
    tspan_
    deltat_::Float64
    p_
    u0_
    rn_
    nutrients_
    eventList_
    envList_ #::Vector{String}()
end

function parseParameters(text)
    pairList = Vector{Pair{Symbol,Float64}}()
    cols = split(text,"\t")
    if(length(cols) < 2)
        return(nothing)
    end
    for col in cols[2:length(cols)]
        col = strip(col)
        columns = split(col,"=")
        if(length(columns) != 2)
            println("ERROR")
            return(nothing)
        end
        symbol = Symbol(strip(columns[1]))
        value = parse(Float64,strip(columns[2]))
        push!(pairList,Pair(symbol,value))
    end
    return(pairList)
end

function parseResources(line)
    text = replace(line,"resources:" => "")
    text = strip(text)
    index1 = findfirst("{",text)
    if(index1 == nothing)
        write(Base.stderr,"\nInvalid input:" * line * "\n\n")
        return(nothing)
    end
    index2 = findfirst("}",text)
    if(index2 == nothing)
        write(Base.stderr,"\nInvalid input:" * line * "\n\n")
        return(nothing)
    end
    text = SubString(text,index1[1]+1,index2[1]-1)
    text = strip(text)
    resourceList = split(text,";")
    map = Dict{String,Any}()
    for text in resourceList
        columns = split(text,"->")
        if(length(columns) != 2)
            write(Base.stderr,"\nInvalid input:" * line * "\n\n")
            return(nothing)
        end
        cols = split(columns[1],"=")
        if(length(cols) != 2)
            write(Base.stderr,"\nInvalid input:" * line * "\n\n")
            return(nothing)
        end
        variable = cols[1]
        value = Base.parse(Float64,cols[2])
        map[variable] = (variable,value,columns[2],line)
    end
    return(map)
end

function parseEnv(line)
    text = strip(line)
    #println("parseEnv:",text)
    m = match(r"ENV\[(.+)\]\s=\s(.+)",text)
    #println("m:",m)
    if(m == nothing)
        return(false)
    end
    if(length(m) != 2)
        return(false)
    end
    return(true)
end

function parseCell(inputFile)
    list = readTextFile(inputFile)
    list = removeBlankLines(list)
    list = removeCommentLines(list,"#")
    index = 0
    nLines = length(list)
    reactionList = nothing
    nCells = nothing
    nSteps = nothing
    p = nothing
    u0 = nothing
    rn = nothing
    tspan = nothing
    dt = nothing
    nutrientMap = nothing
    eventList = Vector{CellEventData}()
    envList = Vector{String}()
    while index < nLines
        index += 1
        line = strip(list[index])
        if(startswith(line,"rn_lamda = "))
            index2 =  findText("end",list,index)
            if(index2 != nothing)
                reactionList = list[index:index2]
                index = index2
                #println(reactionList)
            end
        elseif(startswith(line,"p\t"))
            p = parseParameters(line)
            if(p == nothing)
                message = "\nError, invalid input>" * line * "\n\n"
                abort__(message)
            end
        elseif(startswith(line,"u0\t"))
            u0 = parseParameters(line)
            if(u0 == nothing)
                message = "\nError, invalid input>" * line * "\n\n"
                abort__(message)
            end
        elseif(startswith(line,"nCells="))
            line = replace(line,"nCells=" => "")
            nCells = parse(Int32,line)
        elseif(startswith(line,"nSteps="))
            line = replace(line,"nSteps=" => "")
            nSteps = parse(Int32,line)
        elseif(startswith(line,"timespan="))
            line = replace(line,"timespan=" => "")
            cols = split(line,",")
            t1 = parse(Int32,cols[1])
            t2 = parse(Int32,cols[2])
            tspan = (t1,t2)
        elseif(startswith(line,"dt="))
            line = replace(line,"dt=" => "")
            dt = parse(Float64,line)
        elseif(startswith(line,"resources:"))
            resourceMap = parseResources(line)
            if(resourceMap == nothing)
                exit(0)
            end
        elseif(startswith(line,"cell_event:"))
            event = parseCellEvent(line)
            push!(eventList,event)
        elseif(startswith(line,"nutrients:"))
            nutrientMap = parseResources(line)
            if(nutrientMap == nothing)
                exit(0)
            end
        elseif(startswith(line,"ENV["))
            if(parseEnv(line))
                push!(envList,strip(line))
            end
        else
            message = "\nError, invalid input>" * line * "\n\n"
            abort__(message)
        end
    end
    rn = generateReactionSystem(reactionList)
    #println("RN:",rn)
    params =  model_params(nSteps,nCells,tspan,dt,p,u0,rn,nutrientMap,eventList,envList)
    return(params)
end

function generateReactionSystem(reactionList)
    if(reactionList == nothing)
        return(nothing)
    end
    code = join(reactionList,"\n")
    #println("generateReactionSystem:",code)
    ex = Meta.parse(code)
    eval(ex)
    return(rn_lamda)
end

function xxx_test_event(inputFile)
    println("test_event:",inputFile)
    list = readTextFile(inputFile)
    list = removeBlankLines(list)
    list = removeCommentLines(list,"#")
    nLines = length(list)
    index = 0
    while index < nLines
        index += 1
        line = strip(list[index])
        println(line)
        if(startswith(line,"Events"))
            index2 =  findText("end",list,index)
            println("index:",index2)
            if(index2 != nothing)
                Events.parseEvents(list[index:index2])
                exit()
            end
            break
        end
    end

    println(list)
    return
end


#test_event("event.dat")
