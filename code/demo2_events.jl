

function alwaysUpdate(model,cell,event)
    return(true)
end

function dontSave(model,cell,event)
    return()
end

# demo 2
function calculateGrowth(model,cell,event)
    println("--calculateGrowth")
    #λ = ((c_q(t) + c_m(t) + c_t(t) + c_r(t)) * (gmax*a(t)/(Kgamma + a(t))))/M
    #λ = ((c_q + c_m + c_t + c_r) * (gmax*a/(Kgamma + a)))/M

    eventData = getfield(event,:data_)
    equation = eventData["growthRateEqn"]
    map = Events.getCellParameters(cell)
    growthRate = Events.evaluateEquation(equation,map)
    #result = Events.evaluate__(equation,map)
    growthProg, varType = eventData["growthProg"]
    growthProg += model.dt * growthRate
    eventData["growthProg"] = (growthProg, varType)
    return
end

function calculateGrowthRate(model,cell,event)
    println("calculateGrowthRate")
    eventData = getfield(event,:data_)
    println(eventData)
    exit(999)
end

function checkGrowth(model,cell,event)
    #println("checkGrowth")
    eventData = getfield(event,:data_)
    equation = eventData["growthRateEqn"]
    map = Events.getCellParameters(cell)
    growthLimit, varType = eventData["growthLimit"]
    growthRate = Events.evaluateEquation(equation,map)
    increase = model.dt * growthRate
    growthProg, varType = eventData["growthProg"]
    growthProg += increase
    eventData["growthProg"] = (growthProg, varType)
    growthIncrease, varType = eventData["growthIncrease"]
    eventData["growthIncrease"] = (increase, varType)
    if(growthProg > growthLimit)
        return(true)
    end
    return(false)
end

function resetGrowthEvent(event)
    println("resetGrowthEvent")
    eventData = getfield(event,:data_)
    #println(eventData)
    growthIncrease, varType = eventData["growthIncrease"]
    eventData["growthIncrease"] = (0.0,varType)
    growthProg, varType = eventData["growthProg"]
    eventData["growthProg"] = (0.0,varType)
    return
end

function divideByGrowth(model,cell,event)
    println("divideByGrowth")
    eventData = getfield(event,:data_)
    agent_id = getfield(cell,:id)
    println("DIVIDE CELL:",agent_id, " ", typeof(agent_id))
    nextAgent = model.lastCell + 1
    p_id = getproperty(cell,:id)
    cell_id = getproperty(cell,:id_)
    cellIndex = cell_id - agent_id
    model.lastCell += 2

    #resetGrowthEvent(event)

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

    eventList = getfield(cell1,:events_)
    resetGrowthEvent(eventList[1])
    eventList = getfield(cell2,:events_)
    resetGrowthEvent(eventList[1])
    model.nAgents += 1
println("N Cells:",model.nAgents)
println(length(model.agents))
    #cellSize = sizeof(cell1)
    #println("cell size:",cellSize)
    #cellSize = Base.summarysize(cell1)
    #println("cell size:",cellSize)
    #rn_size = Base.summarysize(cell1.rn_)
    #println("rn size:",rn_size)
    return
end

function hasCellDivided(model,cell,event)
    println("hasCellDivided:",model.nAgents)
    eventData = getfield(event,:data_)
    cellCount, varType = eventData["cellCount"]
    println(cellCount)
    if(model.nAgents > cellCount)
        eventData["cellCount"] = (model,nAgents,varType)
        return(true)
    end
    return(false)
end

function updateNutrient(model,cell,event)
    #println("updateNutrient")
    eventData = getfield(event,:data_)
    nutrient, nutrientType = eventData["nutrient"]
    remainingNutrient, remType = eventData["remainingNutrient"]
    productName, varType = eventData["product"]
    sum = 0.0
    for cell in allagents(model)
        sum += getproperty(cell,productName)
    end
    remainder = nutrient - sum
    eventData["remainingNutrient"] = (remainder,remType)
    return
end
