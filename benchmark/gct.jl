#reader of the GTEx count files, used by the benchmark scripts only (data input stays out of the package)
"""
    read_gct(file)

read a GCT file of read counts (GTEx format, version 1.2 or 1.3)
return a NamedTuple with the counts matrix (genes in rows, samples in columns), the gene names, descriptions and sample ids
"""
function read_gct(file::AbstractString)
    return open(file, "r") do io
        str = readline(io)
        while startswith(str, '#') #dump the version line
            str = readline(io)
        end
        nrows = parse(Int, first(split(str, '\t')))
        header = split(readline(io), '\t')
        offset = header[1] == "id" ? 1 : 0 #version 1.3 has an extra id column before Name and Description
        samples = String.(header[(offset + 3):end])
        counts = zeros(Int, nrows, length(samples))
        genes = fill("", nrows)
        description = fill("", nrows)
        for i in 1:nrows
            fields = split(readline(io), '\t')
            genes[i] = fields[offset + 1]
            description[i] = fields[offset + 2]
            for j in eachindex(samples)
                #some GCT files store counts as floats (e.g. "12.0")
                counts[i, j] = round(Int, parse(Float64, fields[offset + 2 + j]))
            end
        end
        return (counts = counts, genes = genes, description = description, samples = samples)
    end
end
