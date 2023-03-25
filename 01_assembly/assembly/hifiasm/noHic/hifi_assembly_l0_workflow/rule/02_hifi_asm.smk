###################   hifi asm   ######################
rule hifi_asm:
	input:
		fa = "rawdata/02_ccs/{sample}.ccs.fasta.gz",
	output:
		pgfa = "results/02_assembly/02_hifi_asm/{sample}_hifiasm.bp.p_ctg.gfa"
	params:
		dir = "results/02_assembly/02_hifi_asm/{sample}_hifiasm"
	threads: 64
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
			hifiasm -t {threads} -o {params.dir} -l {purge_level} {input.fa}
		"""

rule gfa2fa:
	input:
		primary_gfa = rules.hifi_asm.output
	output:
		primary_fa = "results/02_assembly/02_hifi_asm/{sample}_hifiasm.bp.p_ctg.fa"
	threads: 1
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
			awk '/^S/{{print ">"$2;print $3}}' {input.primary_gfa} > {output.primary_fa}
		"""


rule seqkit_stats:
	input:
		primary_fa = rules.gfa2fa.output.primary_fa,
	output:
		primary_fa_stats = "results/03_assembly_assement/01_seqkit_stats/{sample}_hifiasm.bp.p_ctg",
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
		seqkit stats -a -j 64 -t dna -o {output.primary_fa_stats} {input.primary_fa}
		"""

rule stats_n10_n90:
	input:
		primary_fa = rules.gfa2fa.output.primary_fa,
	output:
		primary_fa_stats = "results/03_assembly_assement/01_seqkit_stats/{sample}_hifiasm.bp.p_ctg_n10_n90",
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
		/home/software/seq_n50.pl {input.primary_fa} > {output.primary_fa_stats}
		"""



rule busco:
	input:
		primary_fa = rules.gfa2fa.output.primary_fa
	output:
		primary_short_summary = "results/03_assembly_assement/02_busco/{sample}_hifiasm.bp.p_ctg/short_summary.specific.solanales_odb10.{sample}_hifiasm.bp.p_ctg.txt"
	shell:
		"""
		cd results/03_assembly_assement/02_busco

			singularity exec ../../../busco.sif busco -m genome \
				-i ../../../{input.primary_fa} \
				-o {Psample}_hifiasm.bp.p_ctg \
				-l ../../../../solanales_odb10 \
				-c 64 -f \
				--offline

		cd ../../../
		"""

