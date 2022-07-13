using CSV, MultipleTesting, DataFrames
cd("/home/pdimens/Omega/USM PhD/Projects/Active/Blackfin Tuna/Analyses/extramiss/fst")
neutral = CSV.read("neutral.hudson.fst", DataFrame)
outlier = CSV.read("outlier.hudson.fst", DataFrame)

function partitionarray(array::AbstractArray, steps::AbstractVector{<:Integer})
    v = axes(array,1)
    v == 1:sum(steps) || error("Steps provided do not sum to length of the first dimension")
    i = firstindex(v)
    tmp = (view(v, i:(i+=s)-1) for s in steps)
    [view(array,r,:) for r in tmp]
end

function adjustpval(fstval::DataFrame)
    fst = deepcopy(fstval)
    rows = size(fst,1)
    pval = mapreduce(vcat, 1:rows-1) do i
        collect(fst[i,i+1:end])
    end
    posthoc = round.(adjust(pval, BenjaminiHochberg()),digits = 4)
    splitpart = partitionarray(posthoc, reverse(collect(1:(rows-1)))) .|> collect
    for i in 1:(rows-1)
        fst[i, (i+1):end] .= splitpart[i][:,1]
    end
    insertcols!(fst, 1, :pop => names(fst))
    return fst
end

CSV.write("neutral.hudson.fdr.fst", adjustpval(neutral))
CSV.write("outlier.hudson.fdr.fst", adjustpval(outlier))
