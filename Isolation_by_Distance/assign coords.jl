using PopGen, DataFrames

data = read_from("bft_maf01_kinless_neutral.gen")
#= 
Locality	Longitude	Latitude
BRZ	-6.438	-34.225
BRZ_SP	0.917214	-29.3466
KEY	24.767	-81.005
MRT	14.586	-61.172
PNS	30.066	-87.214
PR	18.915	-66.536
SCA	32.554	-78.995
TX	27.686	-95.576
VZ	11.491	-66.229
=#

longs = [-6.436, 0.917214, 24.767, 14.586, 30.066, 18.915, 32.554, 27.686, 11.491]
lats = [-34.225, -29.3466, -81.005, -61.172, -87.214, -66.536, -78.995, -95.576, -66.536]

by_pop = groupby(data.meta, :population)

pops = unique(data.meta.population)

for (idx, pop) in enumerate(pops)
    by_pop[(population = pop,)].longitude = longs[idx]
    by_pop[(population = pop,)].latitude = lats[idx]
end
data.meta.latitude = data.meta.latitude |> Vector{Float64}
data.meta.longitude = data.meta.longitude |> Vector{Float64}


# jitter the individual values
for i in 1:length(data.meta.longitude)
    # generate a random number, if it's <0.5 we'll subtract the jitter from the value, otherwise add
    coeff = first(rand(1)) > 0.5 ? 1 : -1
    # same thing as above, except for dividing by 1000 or 10000
    coeff2 = first(rand(1)) > 0.5 ? 0.001 : 0.0001
    data.meta.longitude[i] += coeff * first(rand(1)) * coeff2
    # repeat of the coeff's above but for latitude
    coeff3 = first(rand(1)) > 0.5 ? 1 : -1
    coeff4 = first(rand(1)) > 0.5 ? 0.001 : 0.0001
    data.meta.latitude[i] += coeff3 * first(rand(1)) * coeff4
end

write_to(data, filename = "bft_maf01_kinless_neutral_IDB.gen", format = "ibd")