#
# Interface to NGSlib.R get_annot_tables
#
# prints annotation table to excel
#
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])))
source(file.path(script_dir, "NGSlib.R"))
library(openxlsx);
library(reshape);

usage	= paste("USAGE: Rscript print_annot_table annot_table.tsv annot_table.xlsx [toptaxrank taxgroups]",sep="\n");

args 	= commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
	stop(paste("Missing input files\n",usage,sep="\n"), call.=FALSE);
}

annot_file		= args[1];
excel_file		= args[2];
toptaxrank		= "division";
if(length(args)>=3){
	toptaxrank	= args[3];
}
taxgroups		= c();
if(length(args)>=4){
	taxgroups	= unlist(strsplit(args[4],split = ",",fixed = T));
}

cat(sprintf("\n# \tRscript: get_annot_tables(%s, %s, toptaxrank=%s, taxgroups=c(%s) )\n",
				annot_file,
				excel_file,
				toptaxrank,
				paste(taxgroups,sep=' ',collapse=' ')));

annot		= get_annot_tables(annot_file, excel_file, toptaxrank, taxgroups);
