using Plots, StatsBase, Statistics, LaTeXStrings, JLD, StatsPlots

function FindBorder(dim, Network)
    border0 = []
    for i in 1:dim
        for j in 1:dim
            check = false
            if Network[i,j] == 0
                Neighbors = ((i-1,j),(i,j-1),(i+1,j),(i,j+1))
                for nei in Neighbors
                    if dim+1 ∉ nei && 0 ∉ nei && Network[nei...] == 1
                        check = true
                        break
                    end
                end
                if check == true
                    push!(border0,(i,j))
                end
            end
        end
    end
    return border0
end

function FindCorrelationLength(Network,dim)
    FindOnSpots = findall(x->x==1,Network)
    SpotNums = length(FindOnSpots)
    ΣR² = 0.0; iMC = 0.0; jMC = 0.0
    for indx in FindOnSpots
        iMC += indx[1]/SpotNums
        jMC += indx[2]/SpotNums
    end
    for indx in FindOnSpots
        ΣR² += ((indx[1]-iMC)^2 + (indx[2]-jMC)^2)/SpotNums
    end
    return sqrt(ΣR²)
end

function ClusterGrowth(dim,P)
    Network = zeros(Int,dim,dim)
    # Network[rand(1:dim^2)] = 1
    Network[round(Int, dim/2), round(Int, dim/2)] = 1
    Borders = FindBorder(dim,Network)
    while length(Borders) != 0
        Borders = FindBorder(dim,Network)
        for spot in Borders
            if P > rand()
                Network[spot...] = 1
            else
                Network[spot...] = -1
            end
        end
    end
    return Network
end

function Linear_fit(X, Y)
    A = [hcat(X) reshape(ones(length(X)), length(X), 1)]
    b = reshape(hcat(Y), length(Y), 1)
    line = (A \ b)
    return line
end

# Network = ClusterGrowth(200,0.56)
# heatmap(Network, legend = nothing)

# runnum = 10000
# for p in [0.5,0.55,0.59]
#     totstepdata = []
#     for i in 1:runnum
#         Network = ClusterGrowth(200,p)
#         push!(totstepdata,(FindCorrelationLength(Network,200),length(findall(x-> x==1,Network))))
#         print("\r$i")
#     end
#     println('#')
#     save("../../Data/Q3/Q3-$p-CG.jld", "data", totstepdata)
# end

# AllXiData = []
# AllSData = []
# for p in [0.5,0.55,0.59]
#     Data = load("../../Data/Q3/Q3-$p-CG.jld")["data"]
#     xidata = []
#     sdata = []
#     for data in Data
#         push!(xidata, data[1])
#         push!(sdata, data[2])
#     end
#     push!(AllXiData,xidata)
#     push!(AllSData,sdata)
# end
# AllXiData

#Ploting:
#=
violin(["0.59"],AllXiData[3])
violin!(["0.55"],AllXiData[2])
violin!(["0.50"],AllXiData[1])
plot!(legend = nothing, title = L"Correlation\ Length\ for\ 200\times200\ Network,\ 10000 runs",xlabel = L"P",ylabel=L"\xi_{(P)}")
savefig("../../Figs/Q3/Q3-xi-violon.pdf")


histogram(AllXiData[3],bins = 50, label = L"P=0.59")
histogram!(AllXiData[2],bins = 50, label = L"P=0.55")
histogram!(AllXiData[1],bins = 25, label = L"P=0.50")
plot!(title = L"Correlation\ Length\ for\ 200\times200\ Network,\ 10000 runs",xlabel = L"\xi_{(P)}",ylabel=L"Numbers")
savefig("../../Figs/Q3/Q3-xi-hist.pdf")


violin(["0.59"],AllSData[3])
violin!(["0.55"],AllSData[2])
violin!(["0.50"],AllSData[1])
plot!(legend = nothing, title = L"Cluster\ Size\ for\ 200\times200\ Network,\ 10000 runs",xlabel = L"P",ylabel=L"\S_{(P)}")
savefig("../../Figs/Q3/Q3-S-violon.pdf")
=#

histogram(AllXiData[3],bins = 50, label = L"P=0.59")
histogram!(AllXiData[2],bins = 50, label = L"P=0.55")
histogram!(AllXiData[1],bins = 25, label = L"P=0.50")
plot!(title = L"Cluster\ Size\ for\ 200\times200\ Network,\ 10000 runs",xlabel = L"\S_{(P)}",ylabel=L"Numbers")
savefig("../../Figs/Q3/Q3-S-hist.pdf")


runnum = 10000
AllData = []
for p in 0.5:0.005:0.65
    totstepdata = []
    for i in 1:runnum
        Network = ClusterGrowth(100,p)
        push!(totstepdata,(FindCorrelationLength(Network,100),length(findall(x-> x==1,Network))))
        print("\r$i--$p    ")
    end
    push!(AllData, totstepdata)
end
save("../../Data/Q3/Q3-xi-S-CG.jld", "data", AllData)

AllData
AvgSList

AvgSList = []
STDSList = []
AvgXiList = []
STDXiList = []
for data in AllData
    temspdata = []
    temxipdata = []
    for tup in data
        push!(temxipdata,tup[1])
        push!(temspdata,tup[2])
    end
    push!(AvgSList,mean(temspdata))
    push!(STDSList,std(temspdata))
    push!(AvgXiList,mean(temxipdata))
    push!(STDXiList,std(temxipdata))
end

XData = []
YData = []
for data in AllData
    for tup in data
        push!(XData,tup[2])
        push!(YData,tup[1])
    end
end

scatter(log.(XData), log.(YData),markersize = 2,
    alpha = 0.01, markerstrokewidth=0)
scatter!(title = L"Total\ Data\ (0.5>P>0.65,\ 10000\ runs\ for\ each\ P)",
    legend = nothing, xlabel = L"\xi_{(P)}", ylabel = L"S_{(P)}")
savefig("../../Figs/Q3/Q3-S-XI.pdf")

Line = Linear_fit(log.(AvgSList)[1:end-3],log.(AvgXiList)[1:end-3])
X = hcat(log.(AvgSList)[1]:(log.(AvgSList)[end]-log.(AvgSList)[1])/20:log.(AvgSList)[end])
Y = X .* Line[1,1] .+ Line[2,1]

plot(X,Y, label = L"Y = %$(round(Line[1],digits=3))\dot X + %$(round(Line[2],digits=3))", line = :dash, c = :black)
scatter!(log.(AvgSList), log.(AvgXiList),c = :steelblue, label = L"Data\ Point")
scatter!(title = L"Average\ Data\ (0.5>P>0.65,\ 10000\ runs\ for\ each\ P)",
    xlabel = L"\xi_{(P)}", ylabel = L"S_{(P)}")
savefig("../../Figs/Q3/Q3-S-XI-Avg.pdf")
