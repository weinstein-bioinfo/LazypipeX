#
# INTERFACE TO NGSlib.R QUALITY CONTOR PLOTS
#


usage	= paste("USAGE: Rscript qc_plots.R type file_in file_out",
		"name    [str]  : type of plot, seqlen|readn",
		"file_in [file] : input file, for different plots use different files:",
		"     seqlen : read lengths sep by newlines|fasta file",
		"     readn  : assembly.stats.txt file with reads,reads.flt,reads.rgflt,reads.ass.bwa,reads.genes fields",
		"file_out[file] : file for printing, e.g. readn.jpeg\n",
		sep="\n");

script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])))
source(file.path(script_dir, "NGSlib.R"))


args 	= commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
	stop(paste("missing input\n",usage,sep="\n"), call.=FALSE);
}

type		= tolower(args[1]);
file_in		= args[2];
file_out	= args[3];

if(type == "seqlen"){
	seqlen_hist(file_in = file_in, file_out = file_out);
}
if(type == "readn"){
	name = gsub(".(\\w*)$","",file_in,perl=T);
	plot_readcount(filelist = c(file_in), file_out = file_out, plotpc = T, title=name);
}
