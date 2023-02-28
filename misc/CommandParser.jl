using Pkg

#include("Cell.jl")

module CommandParser

    struct Parameter
        flag_
        type_
        value_
        optionList_
        required_
    end

    function validateParameterValue(parameter,value)
        if (parameter.optionList_ != nothing)
            valid = false
            for i in 2:length(parameter.optionList_)
                println(parameter.optionList_[i])
                if(parameter.optionList_[i] == "*")
                    return(value)
                end
                if( value == parameter.optionList_[i])
                    valid = true
                    break
                end
            end
            println("option list:",parameter.optionList_)
            if(!valid)
                error = "Error, " * parameter.flag_ * ", invalid option:" * value
                CommandParser.__abort(error)
            end
        end
        value = Base.parse(parameter.type_,value)
        return(value)
    end

    function __readFile(inputFile)
        nLines = countlines(inputFile)
        list = String[]
        f = open(inputFile,"r")
        for line in eachline(f)
            if(length(strip(line)) > 0)
                push!(list,String(line))
            end
        end
        close(f)
        return(list)
    end

    function readParameters(inputFile,delimiter="\t")
        unspecified = "*"
        if(!isfile(inputFile))
            error = "Error, invalid file:" * inputFile
            __abort(error)
        end
        list = CommandParser.__readFile(inputFile)
        list = CommandParser.removeCommentLines(list)
        dict = Dict{String,Parameter}()
        for line in list
            parameter = __parseParameter(line,delimiter)
            dict[parameter.flag_] = parameter
        end
        return(dict)
    end

    function __parseParameter(line,delimiter)
        println("__parseParameter:",line)
        cols = split(strip(line),delimiter)
        flag = cols[1]
        type = cols[2]
        value = cols[3]
        default = cols[4]
        required = cols[5]
        println("cols:",cols)
        println(cols[2])
        type = __assignType(cols[2])
        value = __assignValue(type,cols[3])
        optionList =  __assignDefault(type,cols[4])
        required = __assignRequired(cols[5])
        println("OL:",optionList)
        return(Parameter(flag,type,value,optionList,required))
    end

    function __assignType(type)
        type = strip(type)
        if(type == "file")
           return(type)
        elseif(type == "string")
           return(String)
        elseif(type == "int")
           return(Int64)
        else
           println("Error, assignType:",type)
           exit(999)
        end
        return("?")
    end

    function __assignValue(type,value)
        value = strip(value)
        if(value == "none")
            value = "none"
        elseif(value == "*")
            value = nothing
        end
        return(value)
    end

    function __assignDefault(type,value)
        println("assignDefault:",type, " ", value)
        cols = split(strip(value),":")
        println(cols)
        value = cols[1]
        options = split(strip(cols[2]),",")
        if(value == "*")
            return(nothing)
        end
        if(length(cols[2]) == 0)
            return(nothing)
        end
        println(value,"::",options)
        defList = Any[value]
        for option in options
            option = strip(option)
            push!(defList,option)
        end
        #if(type == "file")
        #    println("f")
        #    println(cols)
        #    println(length(cols[2]))
        #    exit()
        #end
        return(defList)
    end

    function __assignRequired(value)
        value = strip(value)
        if(value == "r")
            return(true)
        end
        return(false)
    end

    function __abort(message,errorCode=999)
        write(Base.stderr,"\n" * message * "\n\n")
        exit(errorCode)
    end

    function parse(parameterDictionary,argList)
        println("parse:",argList)
        nArgs = length(argList)
        index = 0
        dict = Dict{String,Any}()
        while(index < nArgs)
            index += 1
            arg = argList[index]
            println(index,":",arg)
            if (!haskey(parameterDictionary,arg))
                error = "Error, invalid arg:" * arg * "\n...aborting"
                CommandParser.__abort(error)
            end
            parameter = parameterDictionary[arg]
            if(parameter.required_)
                index += 1
                if(index > nArgs)
                    error = "Error, " * arg * ", value missing"
                    CommandParser.__abort(error)
                else
                    value = argList[index]
                end
            end
            println("parameter:",parameter)
            value = CommandParser.validateParameterValue(parameter,value)
            println(arg, " Val:",value)
            println("v:",value, " ", typeof(value))
            println(parameter.type_)
            dict[arg] = value
        end
        return(dict)
    end

    function removeCommentLines(textList,commentStr="#")
        list = String[]
        for line in textList
            text = strip(line)
            if(!startswith(text,commentStr))
                push!(list,line)
            end
        end
        return(list)
    end

    function removeComment(text,commentChar)
        index = findfirst(commentChar,text)
        println("index:",index)
        comment = SubString(text,index[1])
        str = SubString(text,1,index[1]-1)
        return(str)
    end

end # module CommandParser

function runCell(dataFile,outputDir)
    ##ENV["GKSwstype"] = "100"
    Random.seed!(6549)
    createDirectory(outputDir)
    model_params = parseCell(dataFile)
    initialiseEnvironment(model_params.envList_)
    #println("env:",ENV)
    model = generateCellModel(model_params)
end

function run(program,dataFile,outputDir,cpus)
    println("run:",program, " ncpu:",cpus) 
    if(program == "cell")
        println("cell")
        runCell(dataFile,outputDir)
    end
end

function test()
    println("test")
    dataFile = "m.dat"
    outputDir = "./model"
    parameterDictionary = CommandParser.readParameters("commandParser.dat")
    println(parameterDictionary)
    dict = CommandParser.parse(parameterDictionary,ARGS)
    println(dict)
    program = nothing
    cpus = 1
    if (haskey(dict,"-run"))
       program = dict["-run"]
    end
    if (haskey(dict,"-cpus"))
       cpus = dict["-cpus"]
    end
    run(program,dataFile,outputDir,cpus)
end

#function main()
#    println("main")
#end

#test()
#main()


