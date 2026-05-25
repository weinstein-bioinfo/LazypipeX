#
# PROJECT:	LAZYPIPE
# SUBJECT:	Interface to histgoram plotting
# DATE: 	2024/02/26

usage	= paste("USAGE: Rscript hist_figures.R type dir_out file_in [file_in2]",
		"type    [str]  : figure type, options: cont_hist/read_hist",
		"dir_out[dir]	: directory for writing histogram figures\n",
		"file_in [file] : sequence length for contig or read1\n",
		"file_in2[file]	: sequence length for read2\n",
		sep="\n");
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])))
source(file.path(script_dir, "NGSlib.R"))

args 	= commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
	stop(paste("missing input\n",usage,sep="\n"), call.=FALSE);
}

type			= tolower(args[1]);
dir_out		= args[2];
file_in		= args[3];
file_in2		= NA;
if(length(args)==4){
	file_in2	= args[4];
}

if(type == "cont_hist"){
	contig_hists(file_in=file_in, outdir=dir_out, outpref= "contig.hist", title="", binwidth=20);
}
if(type == "read_hist"){
	read_hists(read1=file_in, read2=file_in2, outdir=dir_out, outpref="read.hist", binwidth=10);
}
