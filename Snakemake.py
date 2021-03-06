import pdb 

# wildcard_constraints:
#     # sample = "\w{6}-s\w+",
#     group = "\w{6}_\w+"

rule all:
	input: 
		expand("{sample}/miso/summary/miso.miso_summary", sample=config["lookup"].values()),
		config["gff"] + config["feature"] + "/genes_to_filenames.shelve",
		expand(config["compare"] + "{group}/miso_vs_miso/bayes-factors/bygene.miso_bf", group=config["master"]),
		expand(config["compare"] + "{group}/miso_vs_miso/bayes-factors/bygene_b100.0.miso_bf", group=config["master"])



# pdb.set_trace()

rule make_gff_annotation:
	""""Make custom GFF file annotation from refflat file. Script looks for refGene.txt in supplied directory"""
	output:
		config["gff"] + "commonshortest/" + config["feature"] + "." + config["genome_label"] + ".gff3",
		config["gff"] + "commonshortest/" + "gene_lookup." + config["genome_label"] + ".gff3"
	input:
		refgene_dir = "/lab/bartel4_ata/jamie/projects/misoenv/rnaseqlib/",
		gff_dir = config["gff"], 
	params:
		genome = config["genome_label"], 
		other_params = "--flanking-rule commonshortest"
	shell:
		"set +u; source /lab/bartel4_ata/jamie/projects/misoenv/bin/activate; set -u;"
		"python /lab/bartel4_ata/jamie/projects/misoenv/lib/python2.7/site-packages/rnaseqlib-0.1-py2.7.egg/rnaseqlib/gff/gff_make_annotation.py "
		"{input.refgene_dir} {input.gff_dir} --genome-label {params.genome} {params.other_params}"


rule index_gff:
	"""Index GFF file for canonical MISO events"""
	output:
		config["gff"] + config["feature"] + "/genes.gff",
		config["gff"] + config["feature"] + "/compressed_ids_to_genes.shelve", 
		config["gff"] + config["feature"] + "/genes_to_filenames.shelve"
	input:
		gff = config["gff"] + "commonshortest/" + config["feature"] + "." + config["genome_label"] + ".gff3"
	params:
		indexed_gff_dir = config["gff"] + config["feature"]
	shell:
		"set +u; source /lab/bartel4_ata/jamie/projects/misoenv/bin/activate; set -u; "
		"index_gff --index {input.gff} {params.indexed_gff_dir}"


rule index_bam:
	"""Index bam file"""
	output:
		"{sample}/star/Aligned.sortedByCoord.out.bam.bai"
	input:
		"{sample}/star/Aligned.sortedByCoord.out.bam"
	shell:
		"samtools index {input}"


rule run_miso:
	"""Run MISO scripts"""
	output:
		miso_dir = "{sample}/miso",
		has_finished = touch("{sample}/miso/mytask.done")
		# log = "{sample}/miso/logs/main.2018-01-19_{time}.log"
	input:
		is_indexed = config["gff"] + config["feature"] + "/genes.gff", # necessary for dependence on index_gff script
		gff = config["gff"] + config["feature"],
		bam = "{sample}/star/Aligned.sortedByCoord.out.bam",
		bai = "{sample}/star/Aligned.sortedByCoord.out.bam.bai"
	params:
		readlen = "40"
	shell:
		"set +u; source /lab/bartel4_ata/jamie/projects/misoenv/bin/activate; set -u; "
		"miso --run {input.gff} {input.bam} --output-dir {output.miso_dir} --read-len {params.readlen} --use-cluster"


rule summarize_miso:
	"""Summarize MISO output into a single file"""
	output:
		"{sample}/miso/summary/miso.miso_summary",
		has_finished = touch("{sample}/miso/summarize_miso.done")
	input:
		"{sample}/miso/mytask.done" # necessary for dependence on run_miso
	params:
		miso_dir = "{sample}/miso/"
	shell:
		"set +u; source /lab/bartel4_ata/jamie/projects/misoenv/bin/activate; set -u; "
		"summarize_miso --summarize-samples {params.miso_dir} {params.miso_dir}"


def group_samples(wildcards):
	#input: takes in the names of the samples
	#output: needs to return pairs of files with 'eluate_time'
    return [config["lookup"][method + "_" + wildcard] + "/miso/" for method in config["method"] for wildcard in wildcards]

def group_samples_report_file(wildcards):
	#returns the summary file so that compare_miso depends on the out put of summarize_miso
    return [config["lookup"][method + "_" + wildcard] + "/miso/summary/miso.miso_summary" for method in config["method"] for wildcard in wildcards]



rule compare_miso:
	output:
		outdir = config["compare"] + "{group}/miso_vs_miso/bayes-factors/miso_vs_miso.miso_bf",
	input:
		files = group_samples_report_file #for dependence
	params:
		dirs = group_samples # from sample path to group path
	shell:
		"set +u; source /lab/bartel4_ata/jamie/projects/misoenv/bin/activate; set -u; "
		"compare_miso --compare-samples {params.dirs} {output.outfile}"


rule reformat_compare_miso:
	output:
		all_genes = config["compare"] + "{group}/miso_vs_miso/bayes-factors/bygene.miso_bf",
		filtered = config["compare"] + "{group}/miso_vs_miso/bayes-factors/bygene_b100.0.miso_bf"
	input:
		miso_bf = config["compare"] + "{group}/miso_vs_miso/bayes-factors/miso_vs_miso.miso_bf",
		gff_lookup = config["gff"] + "commonshortest/" + "gene_lookup." + config["genome_label"] + ".gff3"
	shell:
		"python /lab/bartel4_ata/jamie/projects/miso_zga/convert_compare_miso.py -g {input.gff_lookup} -c {input.miso_bf} -o {output.all_genes}"

