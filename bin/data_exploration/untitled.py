# https://github.com/kundajelab/chrombpnet/blob/master/chrombpnet/helpers/preprocessing/reads_to_bigwig.py

def bam_to_tagalign_stream(bam_path):
    p = subprocess.Popen(["bedtools", "bamtobed", "-i", bam_path], stdout=subprocess.PIPE)
    return p

def generate_bigwig(
    input_bam_file, 
    input_fragment_file, 
    input_tagalign_file, 
    output_prefix, 
    genome_fasta_file, 
    bsort, 
    tmpdir, 
    no_st, 
    chrom_sizes_file, 
    plus_shift_delta, 
    minus_shift_delta
):
    assert (input_bam_file is None) + (input_fragment_file is None) + (input_tagalign_file is None) == 2, "Only one input file!"

    if input_bam_file:
        p1 = auto_shift_detect.bam_to_tagalign_stream(input_bam_file)
    elif input_fragment_file:
        p1 = auto_shift_detect.fragment_to_tagalign_stream(input_fragment_file)
    elif input_tagalign_file:
        p1 = auto_shift_detect.tagalign_stream(input_tagalign_file)

    if tmpdir is None:
        if bsort:
            cmd = """awk -v OFS="\\t" '{{if ($6=="+"){{print $1,$2{0:+},$3,$4,$5,$6}} else if ($6=="-") {{print $1,$2,$3{1:+},$4,$5,$6}}}}' | sort -k1,1 | bedtools genomecov -bg -5 -i stdin -g {2} | bedtools sort -i stdin """.format(plus_shift_delta, minus_shift_delta, chrom_sizes_file)
        else:
            cmd = """awk -v OFS="\\t" '{{if ($6=="+"){{print $1,$2{0:+},$3,$4,$5,$6}} else if ($6=="-") {{print $1,$2,$3{1:+},$4,$5,$6}}}}' | sort -k1,1 | bedtools genomecov -bg -5 -i stdin -g {2} | LC_COLLATE="C" sort -k1,1 -k2,2n """.format(plus_shift_delta, minus_shift_delta, chrom_sizes_file)
    else:
        assert(os.path.isdir(tmpdir)) # tmp dir path does not exsist
        if bsort:
            cmd = """awk -v OFS="\\t" '{{if ($6=="+"){{print $1,$2{0:+},$3,$4,$5,$6}} else if ($6=="-") {{print $1,$2,$3{1:+},$4,$5,$6}}}}' | sort -T {3} -k1,1 | bedtools genomecov -bg -5 -i stdin -g {2} | bedtools sort -i stdin """.format(plus_shift_delta, minus_shift_delta, chrom_sizes_file, tmpdir)
        else:
            cmd = """awk -v OFS="\\t" '{{if ($6=="+"){{print $1,$2{0:+},$3,$4,$5,$6}} else if ($6=="-") {{print $1,$2,$3{1:+},$4,$5,$6}}}}' | sort -T {3} -k1,1 | bedtools genomecov -bg -5 -i stdin -g {2} | LC_COLLATE="C" sort -T {3} -k1,1 -k2,2n """.format(plus_shift_delta, minus_shift_delta, chrom_sizes_file, tmpdir)

    print(cmd)


    tmp_bedgraph = tempfile.NamedTemporaryFile()
    if no_st:
        print("Making BedGraph (Do not filter chromosomes not in reference fasta)")

        with open(tmp_bedgraph.name, 'w') as f:
            p2 = subprocess.Popen([cmd], stdin=p1.stdout, stdout=f, shell=True)
            p1.stdout.close()
            p2.communicate()
    else:
        print("Making BedGraph (Filter chromosomes not in reference fasta)")

        with open(tmp_bedgraph.name, 'w') as f:
            p2 = subprocess.Popen([cmd], stdin=subprocess.PIPE, stdout=f, shell=True)
            auto_shift_detect.stream_filtered_tagaligns(p1, genome_fasta_file, p2)
            p2.communicate()

    print("Making Bigwig")
    subprocess.run(["bedGraphToBigWig", tmp_bedgraph.name, chrom_sizes_file, output_prefix + "_unstranded.bw"])

    tmp_bedgraph.close()
    
ref_motifs_file=get_default_data_path(DefaultDataFile.atac_ref_motifs)

plus_shift, minus_shift = auto_shift_detect.compute_shift(
    args.input_bam_file,
    args.input_fragment_file,
    args.input_tagalign_file,
    args.num_samples,
    args.genome,
    args.data_type,
    ref_motifs_file
)

print("Current estimated shift: {:+}/{:+}".format(plus_shift, minus_shift))
plus_shift_delta, minus_shift_delta = 4-plus_shift, -4-minus_shift

generate_bigwig(
    args.input_bam_file,
    args.input_fragment_file,
    args.input_tagalign_file,
    args.output_prefix,
    args.genome,
    args.bsort,
    args.tmpdir,
    args.no_st,
    args.chrom_sizes,
    plus_shift_delta,
    minus_shift_delta
)

