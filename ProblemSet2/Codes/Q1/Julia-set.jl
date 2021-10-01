using Plots

cList = [complex(-0.4,-0.6), complex(0,-1),
            complex(-0.12,-0.75), complex(-0.6,0),
                complex(-0.8,0.16), complex(-0.4,0.6)]

function JuliaSet(c)
    abst = zeros((601,601))
    k = (1 + sqrt(1 + (4 * abs(c)))) / 2
    inum = 0
    for i in -1.5:0.005:1.5
        inum += 1
        jnum = 0
        for j in -1.5:0.005:1.5
            jnum += 1
            z = complex(i,j)
            for n in 1:70
                z = (z^2) + c
                if abs(z) > k
                    abst[inum,jnum] = abs(z)
                    break
                end
            end
        end
    end
    return abst
end

for C in cList
    Ju
end
