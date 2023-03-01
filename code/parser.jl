
function abort__(line,message)
    write(Base.stderr,"\n" * line * "\n")
    write(Base.stderr,message)
    write(Base.stderr,"\n")
    exit()
end

function abort__(message)
    write(Base.stderr,"\n")
    write(Base.stderr,message)
    write(Base.stderr,"\n")
    exit()
end

function readTextFile(inputFile)
    if(!isfile(inputFile))
        error = "Error, invalid file:" * inputFile
        println(error)
        exit()
        return(nothing)
    end
    nLines = countlines(inputFile)
    list = String[]
    f = open(inputFile,"r")
    for line in eachline(f)
        push!(list,String(line))
    end
    close(f)
    return(list)
end

function removeBlankLines(textList)
    list = String[]
    for line in textList
        text = strip(line)
        if(length(text) > 0)
            push!(list,line)
        end
    end
    return(list)
end

function removeCommentLines(textList,commentStr)
    list = String[]
    for line in textList
        text = strip(line)
        if(!startswith(text,commentStr))
            push!(list,line)
        end
    end
    return(list)
end

function removeComment(line,commentStr)
    text = strip(line)
    index = findfirst(commentStr,text)
    if(index != nothing)
        text = SubString(text,1,index[1]-1)
        text = rstrip(text)
    end
    return(text)
end

function findText(str,textList,index)
    nLines = length(textList)
    while(index < nLines)
        index +=1
        line = strip(textList[index])
        if(startswith(line,str))
            return(index)
        end
    end
    return(nothing)
end

function parseVariable(line) # variable_name(Float64) = 0.00005
    #println("parseVariable:",line)
    text = strip(line)
    cols = split(strip(text),"=")
    if(length(cols) != 2)
        return(nothing)
    end
    var = strip(cols[1])
    value = strip(cols[2])
    if(!endswith(var,")"))
        return(nothing)
    end
    var = replace(var,")" => "")
    cols = split(var,"(")
    if(length(cols) != 2)
        return(nothing)
    end
    varName = strip(cols[1])
    varType = strip(cols[2])
    dataType = getfield(Main, Symbol(varType))
    #println("varType:",varType, " " , dataType)
    if( startswith(value,"[")  && endswith(value,"]") )
        valueList = value[2:length(value)-1]
        varList = parseVariableList(valueList,dataType)
        if(varList == nothing)
            write(Base.stderr,line * "\n")
            write(Base.stderr,"Error parsing:" * value * "\n")
            exit()
        end
        return(varName,varList,varType)
    elseif(dataType == String)
        return(varName,value,varType)
    elseif(dataType == Symbol)
        return(varName,Symbol(value),varType)
    else
        value = Base.parse(dataType,value)
        return(varName,value,varType)
    end
end

function parseVariableList(text,dataType)
    text = strip(text)
    cols = split(text,",")
    varList = []
    for i in 1:length(cols)
        try
            value = Base.parse(dataType,cols[i])
            push!(varList,value)
        catch e # ArgumentError
            write(Base.stderr,sprint(showerror,e) * "\n")
            return(nothing)
        end
    end
    return(varList)
end

