using Pkg
Pkg.activate(".")
using CSV, DataFrames, StatsBase, Statistics, Dates, Plots

#Indicies of values that are within cutoff (default=3.5) of the IQR of the median
function inlier_indices(d; cutoff = 3.5)
    m,sp = median(d),iqr(d)
    return (m - cutoff*sp) .< d .< (m + cutoff*sp)
end

#Paths and CSV import
latest_GISAID_date = "2023-10-02"
linage_assignment_filepath = "data/2023-10-02_minimal.csv"
export_filepath = "data/2023-10-02_tall_counts.csv"
data = DataFrame(CSV.File(linage_assignment_filepath, types = Dict(:date => String, :deposited_date => String)))
@show countmap(data.qc_status)

#Options
coverage_min_thresh = 0.9
country_count_thresh = 500
lineage_count_thresh = 200
deposition_delay_cutoff = 150 #days

#Fragile - needs updating if new species start showing up. This happens because the GISAID bulk fasta
#download requires scraping info from the names
not_locations = ["Hamster", "cat", "dog", "env", "gorilla", "mink", "monkey", "northern greater galago", "syrian hamster"]

good_quality  = @. data.qc_status == "good"
high_coverage = @. data.coverage > coverage_min_thresh
full_date     = @. occursin(r"\d{4}-\d{2}-\d{2}", data.date)
full_deposited_date     = @. occursin(r"\d{4}-\d{2}-\d{2}", data.deposited_date)
valid_loc     = [!(split(d,"/")[2] in not_locations) for d in data.seqname]

pass = @. good_quality & high_coverage & full_date & full_deposited_date & valid_loc
@show sum(pass), mean(pass)
data_pass = data[pass,:]
data_pass.date = Date.(data_pass.date, dateformat"y-m-d")
data_pass.deposited_date = Date.(data_pass.deposited_date, dateformat"y-m-d")

#Viz the deposition delay...
ENV["GKSwstype"] = "100" #Let the plot run without an X server
mkpath("plots")
pl = histogram([d.value for d in data_pass.deposited_date .- data_pass.date], label = "Frequency", xlabel = "Deposition delay")
plot!([deposition_delay_cutoff,deposition_delay_cutoff],[0,10^5], color = "red", label = "Cutoff", margin = 1Plots.cm)
savefig(pl, "plots/deposition_delay_histogram.svg")

#And retain only wihin threshold
deposited_within_thresh = [d.value <= deposition_delay_cutoff for d in data_pass.deposited_date .- data_pass.date]
data_pass = data_pass[deposited_within_thresh,:]
latest = Date(latest_GISAID_date)

#Checking for other non-location in the second part of the name - list needs manual inspection
weirdlength(n) = length(split(split(n,"|")[1],"/")) > 4
maybenonhuman = data_pass.seqname[weirdlength.(data_pass.seqname)];
f(n) = split(split(n,"|")[1],"/")[2]
println("This should only be locations: ", String.(sort(union(f.(maybenonhuman)))))

lineages = sort(union(data_pass.lineage))
@assert lineages[1] == "BA.2"

deletions = []
#Removing outlier dates
for lin in lineages
    pango_inds = findall(data_pass.lineage .== lin)
    intdays = [d.value for d in (latest .- data_pass.date[pango_inds])]
    outlier_inds = pango_inds[inlier_indices(intdays) .== 0]
    push!(deletions, (lin, length(outlier_inds)))
    deleteat!(data_pass, outlier_inds)
end

country_counts = countmap(data_pass.country)
data_pass = data_pass[[country_counts[c]>country_count_thresh for c in data_pass.country],:]
lineage_counts = countmap(data_pass.lineage)
data_pass = data_pass[[lineage_counts[c]>lineage_count_thresh for c in data_pass.lineage],:]

#Collapsing counts, and exporting
count_df = sort(combine(groupby(data_pass, [:date, :country, :lineage]), nrow => :count))
CSV.write(export_filepath, count_df)

