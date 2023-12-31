#This can all be run from the REPL, in the base package directory
using Pkg
Pkg.activate(".")
using CSV, DataFrames, StatsBase, Statistics, Dates, Flux

#If you're running this on an NVIDIA GPU:
using CUDA
device = gpu
########################################

#If no NVIDIA GPU#######################
#device = cpu
########################################

count_path = "data/2023-10-02_tall_counts.csv"
count_df = DataFrame(CSV.File(count_path))

lineages = sort(union(count_df.lineage))
countries = sort(union(count_df.country))
dates = collect(minimum(count_df.date):maximum(count_df.date))

indmap(v) = Dict(zip(v,1:length(v)))
lineageinds = indmap(lineages)
countryinds = indmap(countries)
dateinds = indmap(dates)

#We'll need numerical "days" for the LF calculations too - We'll make day zero the most recent date
days = [d.value for d in dates .- maximum(dates)]

#lineages-by-countries-by-days
count_tensor = zeros(Int64,length(lineages), length(countries), length(dates))
d = count_df
for r in 1:size(d,1)
    count_tensor[lineageinds[d.lineage[r]],countryinds[d.country[r]],dateinds[d.date[r]]] = d.count[r]
end

#Setting up model using Flux
struct MLR
  R
  C
end
Flux.@functor MLR
MLR(T, lineages::Integer, countries::Integer) = MLR(zeros(T,lineages), zeros(T,lineages,countries,1))
#Scaling of R makes rates interpretable in years, scaling of C to help model convergence
(m::MLR)(days) = Flux.unsqueeze((m.R ./ 365.0f0)*days', dims = 2) .+ m.C .* 10.0f0

#Init model, and get it onto the GPU
model = MLR(Float32,length(lineages),length(countries)) |> device
days, count_tensor = (days, count_tensor) |> device

#Fitting model
opt = Flux.Optimiser(Flux.WeightDecay(1.0f-10), AdamW(0.1f0))
ps = Flux.params(model)

current_l = Inf
@time for i in 1:100000 #323.014108 seconds
    l,gs = Flux.withgradient(ps) do
            -mean(logsoftmax(model(days)) .* count_tensor)
        end
    Flux.Optimise.update!(opt, ps, gs)
    #Everything below here is monitoring the fit, adjusting the learning rate, and stopping
    if i % 100 == 0
        model.R .= model.R .- mean(model.R) #Centering R, removing an identifiability issue
        print(l, "; ", opt[2][1].eta, "; ")
        #Monitoring growth of 3 lineages spread over the expected range of R
        println([model.R[lineages .== "BA.2"][1], model.R[lineages .== "XBB.1.5"][1], model.R[lineages .== "HK.3"][1]])
        if l > current_l
            opt[2][1].eta = opt[2][1].eta * 0.85f0
            println("LR decay")
            if opt[2][1].eta < 0.0001f0 #Stopping criterion
                println("Done")
                break
            end
        else
            current_l = l
        end
    end
end
#Should look something like:
#0.1746748; 0.00010854139479114348; Float32[-71.824394, 27.58488, 85.06291]

#Saving model fit
using BSON
mkpath("model_fits")
cpumodel = cpu(model) #BSON doesn't like GPU arrays
@BSON.save "model_fits/modelfit.bson" cpumodel

#Saving inferred rates
rate_df = DataFrame()
rate_df.pango = lineages
rate_df.R = Float64.(model.R)
rate_df.seq_volume = sum(count_tensor,dims = (2,3))[:]
CSV.write("model_fits/rates.csv",rate_df)
for c in 1:length(countries)
    rate_df[!,"seq_volume_$(countries[c])"] = sum(count_tensor[:,c,:],dims = 2)[:]
end
CSV.write("model_fits/rates_with_country_volumes.csv",rate_df)

#Plotting fit to countries
using Plots
plotpath = "plots"
mkpath(plotpath)
model, days, count_tensor = (model, days, count_tensor) |> cpu
mod_fit = softmax(model(days), dims = 1)
freq_tensor = count_tensor ./ sum(count_tensor, dims = 1) #This is meant to have a lot of NaNs - don't panic

ENV["GKSwstype"] = "100" #Let the plot run without an X server

for country_to_plot in 1:length(countries)
    topN2scatter = 25
    topNsequenced = sortperm(sum(count_tensor[:,country_to_plot,:], dims = 2)[:], rev = true)[1:topN2scatter]
    point_scal = 100/sqrt(sum(count_tensor[:,country_to_plot,:]))
    pl = plot(size = (900,300))
    color_scale = cgrad(:rainbow, length(lineages), categorical = true);
    for l in 1:length(lineages)
        plot!(days, mod_fit[l,country_to_plot,:], legend = :none, label = :none, color = color_scale[l])
        if l in topNsequenced
            notmissing = .!isnan.(freq_tensor[l,country_to_plot,:][:])
            scatter!(   days[notmissing], 
                        freq_tensor[l,country_to_plot,:][notmissing], 
                        markersize = point_scal .* sqrt.(count_tensor[l,country_to_plot,:][notmissing]),
                        label = lineages[l], markerstrokewidth = 0, legend = :outertopleft,
                        alpha = 0.5, color = color_scale[l])
        end

    end
    savefig(pl, joinpath(plotpath,"$(countries[country_to_plot]).pdf"))
end
