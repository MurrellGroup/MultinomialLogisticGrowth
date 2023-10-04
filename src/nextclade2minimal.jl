#Nextclade command generating the input for this:
#nextclade run --silent --input-dataset BA2_nextclade_refdata --output-csv 2023_10_02.csv GISAID/sequences.fasta

using Pkg
Pkg.activate(".")

using CSV, DataFrames, Dates

function parse_name(seqname)
    virus, date1, date2 = split(seqname, '|')
    names = split(virus, '/')
    country = names[2] == "env" ? names[3] : names[2]
    return seqname, country, date1, date2
end

nextclade_csv_path = "data/2023_10_02.csv"
minimal_export_path = "data/2023-10-02_minimal.csv"
NCdf = CSV.read(nextclade_csv_path, DataFrame, delim = ';')

parsed_names = parse_name.(NCdf.seqName)
@assert [p[1] for p in parsed_names] == NCdf.seqName

#"seqname,country,date,lineage,qc_status,coverage"
minimal_df = DataFrame()
minimal_df.seqname = NCdf.seqName
minimal_df.country = [p[2] for p in parsed_names]
minimal_df.date = [p[3] for p in parsed_names]
minimal_df.deposited_date = [p[4] for p in parsed_names]
minimal_df.lineage = NCdf.Nextclade_pango
minimal_df.qc_status = NCdf."qc.overallStatus"
minimal_df.coverage = NCdf.coverage
minimal_df = minimal_df[.!ismissing.(minimal_df.lineage),:]

CSV.write(minimal_export_path,minimal_df)
