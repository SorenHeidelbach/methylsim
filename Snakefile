

configfile: "config.yml"

rule all:
    input:
        expand("results/ZymoFecal_contig_1330/{snp_rate}/nanomotif_motif_discovery", snp_rate=config["reference_mutation"]["snp_rate"]),


rule mutate_reference:
    input:
        lambda wildcards: config["samples"][wildcards.sample].get("ref", "resources/reference.fasta")
    output: 
        outdir=directory("results/{sample}/{mutation_params}/mutation_simulator"),
        ref="results/{sample}/{mutation_params}/mutation_simulator/mutated_reference_ms.fasta"
    shell:
        """
        mkdir -p {output.outdir}
        pixi run -e mutation-simulator -- \
            mutation-simulator {input} -o {output.outdir}/mutated_reference args --snp {wildcards.mutation_params}
        """


rule simulate_reads_from_mutated_reference:
    input:
        ref="results/{sample}/{mutation_params}/mutation_simulator"
    output:
        reads1="results/{sample}/{mutation_params}/simulated_reads.fastq"
    params:
        motifs=lambda wildcards: config["methylation_sim"].get("motifs", []),
        motif_high_prob_rate=lambda wildcards: config["methylation_sim"]["motif_high_prob_rate"],
        non_motif_high_prob_rate=lambda wildcards: config["methylation_sim"]["non_motif_high_prob_rate"],
        coverage=lambda wildcards: config["methylation_sim"]["coverage"]
    shell:
        """
        pixi run -e methylsim -- \
        methylsim --reference {input.ref}/mutated_reference_ms.fasta --motif {params.motifs} \
         --motif-high-prob {params.motif_high_prob_rate} --non-motif-high-prob {params.non_motif_high_prob_rate} --output-fastq {output} \
         --simulator badreads \
         --badreads-extra  '--quantity {params.coverage}x'
        """


rule map:
    input:
        reads="results/{sample}/{mutation_params}/simulated_reads.fastq",
        ref="results/{sample}/{mutation_params}/mutation_simulator/mutated_reference_ms.fasta"
    output:
        bam="results/{sample}/{mutation_params}/mapped_reads.bam"
    shell:
        """
        minimap2 -ay -x map-ont {input.ref} {input.reads} | samtools sort -o {output.bam}
        samtools index {output.bam}
        """



rule modkit_pileup:
    input:
        bam="results/{sample}/{mutation_params}/mapped_reads.bam"
    output:
        pileup="results/{sample}/{mutation_params}/modkit_pileup.bed"
    shell:
        """
        pixi run -e modkit -- \
        modkit pileup {input.bam} {output.pileup}
        """



rule nanomotif_motif_disocevery:
    input:
        pileup="results/{sample}/{mutation_params}/modkit_pileup.bed",
        ref="results/{sample}/{mutation_params}/mutation_simulator/mutated_reference_ms.fasta",
        contig_bin=lambda wildcards: config["samples"][wildcards.sample].get("contig_bin", "resources/contig_bins.tsv")
    output:
        directory("results/{sample}/{mutation_params}/nanomotif_motif_discovery")
    shell:
        """
        pixi run -e nanomotif -- \
        nanomotif motif_discovery {input.ref} {input.pileup}  -c {input.contig_bin} --out {output} -v
        """


