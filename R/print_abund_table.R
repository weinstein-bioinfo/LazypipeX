#
# Interface to NGSlib.R get_abund_tables
#
# prints abundance table to excel
#
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])))
source(file.path(script_dir, "NGSlib.R"))
library(openxlsx);
library(reshape);

usage	= paste("USAGE: Rscript print_abund_table abund_table.txt abund.xlsx taxranks [taxgroups tail]",sep="\n");

args 	= commandArgs(trailingOnly=TRUE)
if( length(args)<3 ) {
	stop(paste("Missing input files\n",usage,sep="\n"), call.=FALSE);
}

abund_table		= args[1];
excel_file		= args[2];
taxranks			= unlist(strsplit(args[3],split = ",",fixed = T));
taxgroups		= c();
if(length(args)>=4){
	taxgroups	= unlist(strsplit(args[4],split = ",",fixed = T));
}
tail 			= 0;
if(length(args)>=5){
	tail		= as.numeric(args[5]);
}

cat(sprintf("\n# \tRscript: get_abund_tables(%s, %s, taxranks=c(%s), taxgroups=c(%s), tail=%d)\n",
				abund_table,
				excel_file,
				paste(taxranks,sep=',',collapse=','),
				paste(taxgroups,sep=',',collapse=','),
				tail));

abund		= get_abund_tables(abund_table=abund_table, excel_file=excel_file, taxranks=taxranks, taxgroups=taxgroups, tail=tail);
