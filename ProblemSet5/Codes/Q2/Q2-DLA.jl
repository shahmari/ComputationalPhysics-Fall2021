L = 200

#Diffusion-limited aggregation


H1 = 10
H2 = H1 + 5
choices = ((1,0),(0,1),(-1,0),(0,-1))
Network = zeros(L,L,2)
Network[:,1,1] = ones(L)
Network


RWpos = [rand(1:L), H1]

ParticlesNumber = 0
while H2 < L
    RWpos .+= rand(choices)
    if RWpos[2] <= H1 - 5
        for neighbor in ((1,0),(0,1),(-1,0),(0,-1))
            if RWpos[1] + neighbor[1] == 0 && Network[L,RWpos[2] + neighbor[2],1] != 0
                Network[RWpos...,2] += 1
                if Network[RWpos...,2] == 3
                    Network[RWpos...,1] = 1
                end
                break
            elseif RWpos[1] + neighbor[1] == L + 1 && Network[1,RWpos[2] + neighbor[2],1] != 0
    end
end

heatmap(rotr90(Network))
