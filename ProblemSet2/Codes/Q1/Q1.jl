using Plots

cList = [complex(-0.4,-0.6), complex(0,-1),
            complex(-0.12,-0.75), complex(-0.6,0),
                complex(-0.8,0.16), complex(-0.4,0.6)]

function JuliaSet(c, opnum)
    abst = zeros((601,601))
    k = (1 + sqrt(1 + (4 * abs(c)))) / 2
    inum = 0
    for i in -1.5:0.005:1.5
        inum += 1
        jnum = 0
        for j in -1.5:0.005:1.5
            jnum += 1
            z = complex(i,j)
            for n in 1:opnum
                z = (z^2) + c
                if abs(z) > k
                    abst[inum,jnum] = abs(z)
                    break
                # else
                #     abst[inum,jnum] = k
                end
            end
        end
    end
    return abst
end

# theme(:dark)

for C in 1:length(cList)
    data = JuliaSet(cList[C], 100)
    heatmap(data, c = cgrad(:solar, rev = false), legend = false, border=:none, title = "C = $(cList[C])")
    savefig("C:\\Users\\Yaghoub\\Documents\\GitHub\\ComputationalPhysics-Fall2021\\ProblemSet2\\Figs\\Q1\\JuliaSet$C.png")
end
