using PopGen
cd("~/Blackfin Tuna/Analyses/Isolation by Distance/")
infile = "../BFT_nohaplo_neutral.gen"
infile_haplo ="BFT_new_maxmiss80.gen"

#infile_win = "C:/Users/pdime/Desktop/BFT_maxmiss80.txt"
bft = genepop(infile) # numpops removed in updated version of PopGen
#bft = genepop(infile_haplo)
## BRZ KEY LA MRT PNS PR SCA TX VZ
popcounts = population(bft).count
current_popnames = unique(bft.samples.population)
popnames = ["BRZ", "BRZ_SP","KEY", "MRT", "PNS", "PR", "SCA", "TX", "VZ"]
rn_dict = Dict()
[rn_dict[i] = j for (i,j) in zip(current_popnames, popnames)]
population!(bft, rename = rn_dict)

longitudes = [-6.438, 0.917214, 24.767, 14.586, 30.066, 18.915, 32.554, 27.686, 11.491]
latitudes = [-34.225, -29.3466, -81.005, -61.172, -87.214, -66.536, -78.995, -95.576, -66.229]

open("coords.txt", "w") do file
    println(file, "Locality\tLongitude\tLatitude")
    for (i,j,k) in zip(popnames, longitudes, latitudes)
        println(file, i,"\t", j, "\t", k)
    end
end

longitudes_jitt = []
latitudes_jitt = []

for (long, lat, pcount) in zip(longitudes,latitudes, popcounts)
    long_jitt = round.(long .+ (rand(-1.0:0.001:1.0, pcount) .* 0.01), digits = 5)
    latt_jitt = round.(lat .+ (rand(-1.0:0.001:1.0, pcount) .* 0.01), digits = 5)
    push!(longitudes_jitt, long_jitt)
    push!(latitudes_jitt, latt_jitt)
end

longitudes_jitt_flat = Iterators.flatten(longitudes_jitt) |> collect
latitudes_jitt_flat = Iterators.flatten(latitudes_jitt) |> collect
oldnames = bft.samples.name
newnames = Vector{String}()

for (i,j,k) in zip(longitudes_jitt_flat,latitudes_jitt_flat,oldnames)
    push!(newnames, string("$i\t$j\t$k"))
end

# write to output file
open("bft_jitt_loc.txt", "w") do f
    for i in newnames
        println(f, i)
    end
end

popgenfile = open(readlines, infile_haplo)
popgenfile = open(readlines, "C:/Users/pdime/Desktop/BFT_maxmiss80.txt")
outfile ="BF_haplo_IDB.txt"
outfile = "C:/Users/pdime/Desktop/BFT_IBD.txt"
open(outfile, "w") do f
    name_num = 1
    len = length(popgenfile)
    for line in 1:len
        if occursin("POP", popgenfile[line]) == true
        # A little condition to add the first POP to the file
            if name_num == 1
                println(f, "POP")
            end
            # skip POP lines
            continue
        elseif occursin(",",  popgenfile[line]) != true
            # output locus names
            println(f,  popgenfile[line])
            continue
        else
            eachline_mod = replace(popgenfile[line], oldnames[name_num] => newnames[name_num])
            if line != len
                println(f, eachline_mod, "\nPOP")
                name_num += 1
            else
                print(f, "$eachline_mod")
            end
        end    end
end

open(outfile, "w") do f
    name_num = 1
    len = length(popgenfile)
    for line in 1:len
        if line == 1
            println(f, popgenfile[line])
            continue
        elseif line == 2
            println(f, popgenfile[line])
            continue
        end
        if occursin("POP", popgenfile[line]) == true
        # A little condition to add the first POP to the file
            if name_num == 1
                println(f, "POP")
            end
            # skip POP lines
            continue
        elseif occursin(",",  popgenfile[line]) != true
            # output locus names
            println(f,  popgenfile[line])
            continue
        else
            eachline_mod = replace(popgenfile[line], oldnames[name_num] => newnames[name_num])
            if line != len
                println(f, eachline_mod, "\nPOP")
                name_num += 1
            else
                print(f, "$eachline_mod")
            end
        end
    end
end
