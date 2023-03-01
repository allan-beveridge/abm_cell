
function agentPlot(df,nCells,outputDir,outputPrefix)
    plot1_map = Dict([("a","red"),("r","blue")])
    plot2_map = Dict([("m_r","red"),("m_q","blue"),("e_t","green"),("e_m","yellow"),("q","pink"),("m_t","brown"),("m_m","grey"),("c_r","purple"),("c_t","orange"),("c_m","cyan"),("c_q","magenta"),("s_i","black")])
    plot_map_all = Dict([("r","red"),("e_t","blue"),("e_m","green"),("q","yellow"),("m_r","pink"),("m_t","brown"),("m_m","grey"),("m_q","purple"),("c_r","orange"),("c_t","cyan"),("c_m","magenta"),("c_q","black"),("a","crimson"),("s_i","indigo")])
    for i in 1:nCells
        filePrefix = outputPrefix  * "_Cell" * string(i)
        outputFile = outputDir * "/" * filePrefix * ".csv"
        filtered_df = filter(row -> row.id == i, df) # filter main output file by cell id
        println(outputFile)
        CSV.write(outputFile,filtered_df)            # output results for each cell id to separate file
        plotFile = outputDir * "/" * outputPrefix * "_a_r.png"
        plot_results(filtered_df,plot1_map,plotFile)
        plotFile = outputDir * "/" * outputPrefix * "_2.png"
        plot_results(filtered_df,plot2_map,plotFile)
        plotFile = outputDir * "/" * outputPrefix * "_all.png"
        plot_results(filtered_df,plot_map_all,plotFile)
    end
    return
end

function plot_results(df,plotMap,plotFile)
    #display = ENV["DISPLAY"]
    plot(xlabel="Time Steps",ylabel="Concentration")
    for (column,colour) in plotMap
        name = Symbol(column)
        plot!(df.step,df[:,name],color=colour,label=column)
        #    #plot(df.step,df[:,name],xlabel="Time",ylabel="Concentration",color=colour,label=column)
    end
    savefig(plotFile)
    return
end

