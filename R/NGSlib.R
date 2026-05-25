# 
# FUNCTIONS FOR NGS PIPE
#

#
# Sorts files in a list using one of defined rules:
# 'numsuffix'    by numeric suffix, e.g. 'file_B_1' > 'file_A_2' > 'file_A3'
sort.files<- function(filelist, rule='numsuffix'){
	i<- order(as.numeric(gsub("(.*[A-z_-]+)([0-9]+)([A-z]?)$", "\\2", filelist)));
	return(filelist[i]);
}

# Converts filenames to short labels
# by default uses string prefix + the first occuring number
file2label<- function(filelist){
}

# converts data frame to flat format 
# (columnwise catenation with addition of rowname and colname columns)
convert_to_flat<- function(data){
	m<- nrow(data);
	n<- ncol(data);
	data_flat<- as.data.frame(matrix(0, nrow= m*n, ncol=3));
	
	colnames(data_flat)<- c("rowname","colname","value");
	data_rownames<- rownames(data);
	data_colnames<- colnames(data);
	
	for(i in 1:m){
		for(j in 1:n){
			data_flat[(j-1)*m+i,1]<- data_rownames[i];
			data_flat[(j-1)*m+i,2]<- data_colnames[j];
			data_flat[(j-1)*m+i,3]<- data[i,j];
		}
	}
	data_flat$rowname<- factor(data_flat$rowname, levels = unique(data_flat$rowname));
	data_flat$colname<- factor(data_flat$colname, levels = unique(data_flat$colname))
	return(data_flat)
}

# QUALITY CONTROL: READ LENGTH HISTOGRAM
#
seqlen_hist<- function(file_in, file_out = "none", name = "none", width=5.2, res=300, pointsize=9){
	
	# DEFAULTS
	if(file_out == "none"){ file_out<- gsub(".(\\w*)$",".hist1.jpg",file_in,perl=T);}
	if(name == "none"){ name = gsub(".(\\w*)$","",file_in,perl=T);}
	
	if(file_out!="none"){
		jpeg(filename = file_out, width=width*res, height=width*res, units="px", res=300, pointsize=pointsize);
	}     
	# HISTOGRAM FROM LENGTH VECTOR FILE
	if(length(grep(".(len|length)$",file_in,perl=T,ignore.case=T)) > 0){
		readlen<- read.csv(file=file_in, sep="", header=FALSE)
		readlen<- readlen[,1]
		hist(readlen, breaks=seq(0,max(readlen)+10,by=10), main=name, xlab="seq length");
		#par(new=T);
		#plot(density(readlen), main="", ylab="",xlab="",col="red", yaxt="n")
		axis(4);
		
		# mean+median+quantiles
		q<- quantile(readlen,probs = seq(0,1,0.25))
		s<- sprintf("mean\t: %1.1f\nmedian\t: %1.1f\nQ0\t: %1.1f\nQ25\t: %1.1f\nQ50\t: %1.1f\nQ75\t: %1.1f\nQ100\t: %1.1f",mean(readlen),median(readlen),q[1],q[2],q[3],q[4],q[5]);
		mtext(s,side=3,adj = 0.1,padj = 1.1);
	}
	if(file_out!="none"){
		cat("# histogram printed to ",file_out,"\n")
		dev.off();
	}    
}


# Plot read numbers from raw to reads in genes
# filelist      list of assembly.stats.txt files
plot_readcount<- function(filelist, file_out="none",title="",legend=F,labels=basename(dirname(filelist)),
		plotpc = F, width=5.2, res=300, pointsize=9){
	if(!interactive()) pdf(NULL);
	
	if(require(colorspace, quietly=T) == TRUE){
		pal<- rainbow_hcl(n=length(filelist));
	}
	else{
		pal<- rep("#000000",length(filelist));
	}
	
	
	data<- read.table(file=filelist[1],header=F,sep="\t",row.names=1);
	if(length(filelist) >1){
		for(file in filelist[2:length(filelist)]){
			column<- read.table(file = file,header = F,sep="\t",row.names=1);
			data<- cbind(data,column);
		}
	}
	data<- t(data);
	rownames(data)<- 1:nrow(data);
	data<- as.data.frame(data);
	
	# converting survival rates to pcs
	if(plotpc){
		data$reads.flt    <- (data$reads.flt)/data$reads;
		data$reads.hgflt  <- data$reads.hgflt/data$reads;
		data$reads.contigs <- data$reads.contigs/data$reads;
		data$reads.orfs  <- data$reads.orfs/data$reads;
		data$reads        <- data$reads/data$reads;
		
		yticks<- seq(0, 1,by = 0.1);
		ylabels<- sprintf("%.0f",yticks*100);
		ylim<- c(0,1);
		ylab<- "Number of reads (%)";
	}
	else{
		#ysteP <- ceiling(max(data$reads)/10);
		yticks<- seq(0, max(data$reads)*1.05,length.out=11);
		ylabels<- sprintf("%iK",yticks/1000);
		ylim<- c(0, max(data$reads)*1.05);
		ylab<- "Number of reads";
	}
	
	par(mar= c(5,6,4,9)+0.1);
	par(mgp= c(4,1,0), xpd=T);
	
	if(file_out!="none"){
		jpeg(filename = file_out, width=width*res, height=width*res, units="px", res=300, pointsize=pointsize);
	}    
	for(i in 1:nrow(data)){
		if(i>1){ par(new=T)};
		
		plot(c(1:5),data[i,c("reads","reads.flt","reads.hgflt","reads.contigs","reads.orfs")],
				type='b',lty=1,lwd=2,col=pal[i],pch=i,ylab =ylab, xlab="", yaxt="n",xaxt="n",ylim=ylim);
	}
	axis(1, at=c(1:5),labels=c("reads","reads.flt","reads.hgflt","reads.contigs","reads.orfs"), tick=T,las=1)
	axis(2, at=yticks,labels=ylabels, tick=T,las=2)
	if(title != ""){ title(main=title); }
	#if(subtitle != ""){ title(sub=subtitle);}
	if(legend){
		legend("right",legend=labels,lty=1,lwd=2,pch=1:nrow(data),col=pal, bty="n",cex=1,inset=-0.4);
	}
	
	if(file_out!="none"){
		dev.off();
	}
	par(mar= c(5,4,4,2)+0.1, xpd=F);
}

# input:
# stats     assembly stats, including stats$reads
# timer     timer stats, including timer$cpu
plot_timer1<- function(stats,timer){
	ystep = 100000;
	yticks<- seq(0, ceiling(max(stats$reads)/ystep)*ystep, by=ystep);
	ylabels<- sprintf("%iK",yticks/1000);
	ylim<- c(0, ceiling(max(stats$reads)/ystep)*ystep);    
	
	par(mar= c(5,6,4,2)+0.1);
	
	mod1<- lm(stats$reads ~ timer$cpu)
	plot(timer$cpu, stats$reads, type='p',lty=1,lwd=2,col="blue",pch=15,
			xlab="CPU(s)",ylab="Number of reads",yaxt="n",ylim=ylim)
	axis(2, at=yticks,labels=ylabels, tick=T,las=2)
	title(main="CPU performance")
	abline(mod1, lty=1,lwd=2,col="black")
	
	s<- sprintf("reads/s   : %1.1f\n",mod1$coefficients[[2]]);
	mtext(s,side=3,adj = 0.1,padj = 1.5,cex=1.2);
	
	#pre <- predict(mod1) # plot distances between points and the regression line
	#segments(stats$reads,timer$cpu,stats$reads,pre,col="red")
	
	#mod2<- lm(stats$reads ~ timer$real);
	par(mar= c(5,4,4,2)+0.1, xpd=F);
}

# boxplot for time perf across pipeline parts
plot_timer2<- function(timer2,lwd=1){
	
	tmp<- convert_to_flat(timer2);
	colnames(tmp)<- c("sample","pipeline_step","cpu");
	plot(tmp$pipeline_step,tmp$cpu, main="Stepwise CPU",lwd=lwd, xlab="pipeline step", ylab="CPU(%)");
	labels<- c("1:histogramm","2:Trimmomatic filter","3:refgen filter","4:assemble","5:MetaGene",
			"6:SANSparallel","8:realign reads","9:sort by taxa","10:realign within taxa",
			"11:stats","12:collect results","13:cleanup");
	legend("topright",legend=labels,col="black",cex=1,inset=0.05);
}

# EVALUATES PERFORMANCE ON A GIVEN BENCHMARK
# pred_file     : taxa predictions as data.table with fields: <readn,species,species_id,genus,genus_id,family,family_id,superkingdom, superkingdom_id>
# true_file     : true abundance as data.table (same format)
# excel_file    : print resutls to this file
# pred_phageflt : filter phages from pred_file using NGS.phage.filter.R
# true_phageflt : filter phages from true_file using NGS.phage.filter.R
# host_taxid    : reference genome taxid to exclude from predictions and read dist extimation
# tail          : the least abundant taxa that form up to tail% of read distribution will be ignored
#
bmark_stats<- function(pred_file,true_file,excel_file="none",filter_phages_true=F,filter_phages_pred,host_taxid=9606, tail_vi=1, tail_ba=5){
	library(reshape);
	library(openxlsx);
	res<- list();
	
	# READ DATA
	pred<- read.table(pred_file,sep="\t",header=T, comment.char="",stringsAsFactors = F,quote="");
	true<- read.table(true_file, sep="\t",header=T, comment.char = "",stringsAsFactors = F,quote="");
	
	# SORT BY READN
	pred<- pred[order(pred$readn,decreasing=T),];
	true<- true[order(true$readn,decreasing=T),];
	
	# EXCLUDE HOST TAXID
	ind<- which(pred$taxid == host_taxid);
	if(length(ind)>0){
		pred<- pred[-ind,];
	}
	
	# FILTER PHAGES FROM GROUND TRUTH
	if(filter_phages_true){
		ind<- which(true$class %in% phage.classes);
		ind<- union(ind, which(true$order %in% phage.orders));
		ind<- union(ind, which(true$family %in% phage.families));
		ind<- union(ind,grep("\\<phage\\>|\\<bacteriophage\\>",true$species,ignore.case = T));
		true[ind,"superkingdom"]<- "Bacteriophage";
	}
	
	# FILTER PHAGES FROM PREDICTIONS
	if(filter_phages_pred){
		#source("R/NGS.phage.filter.R");
		ind<- which(pred$class %in% phage.classes);
		ind<- union(ind, which(pred$order %in% phage.orders));
		ind<- union(ind, which(pred$family %in% phage.families));
		ind<- union(ind,grep("\\<phage\\>|\\<bacteriophage\\>",pred$species,ignore.case = T));
		pred[ind,"superkingdom"]<- "Bacteriophage";
	}
	
	
	# DELETE EXTRAC FIELDS: THESE WILL INTERFERE WITH MELT()>RECAST IF NOT NUMERICAL 
	pred<- pred[, c("readn","species","species_id","genus","genus_id","family","family_id","superkingdom","superkingdom_id")];
	true<- true[, c("readn","species","species_id","genus","genus_id","family","family_id","superkingdom","superkingdom_id")];
	
	tmp<- melt(pred, id=c("species","species_id","genus","genus_id","family","family_id","superkingdom","superkingdom_id"))
	tmp2<- cast(tmp,species+species_id+genus+genus_id+family+family_id+superkingdom+superkingdom_id~variable, sum)
	pred<- tmp2[order(tmp2$readn,decreasing=T),]    
	
	tmp<- melt(true, id=c("species","species_id","genus","genus_id","family","family_id","superkingdom","superkingdom_id"))
	tmp2<- cast(tmp,species+species_id+genus+genus_id+family+family_id+superkingdom+superkingdom_id~variable, sum)
	true<- tmp2[order(tmp2$readn,decreasing=T),]
	
	# ADD CSUM COL + exclude CSUM TAIL
	readn_sum<- sum(pred$readn);
	pred$csum<- (cumsum(pred$readn)-pred$readn)/readn_sum; # cumsum of preciding taxa
	ind.vi<- which(pred$csum < (100.0-tail_vi)/100.0);
	ind.ba<- which(pred$csum < (100.0-tail_ba)/100.0);
	pred.vi<- pred[ind.vi,];
	pred.ba<- pred[ind.ba,];
	
	
	# saving true, pred, TP, FP and FN for Viruses and Bacteria
	res$vi$true<-   true[true$superkingdom=="Viruses",c("species","species_id","readn")];
	res$vi$pred<-   pred.vi[pred.vi$superkingdom=="Viruses",c("species","species_id","readn")];
	res$vi$TP<-     res$vi$pred[
			match( intersect(res$vi$pred$species_id,res$vi$true$species_id), res$vi$pred$species_id), ];
	res$vi$FP<-     res$vi$pred[
			match( setdiff(res$vi$pred$species_id, res$vi$true$species_id), res$vi$pred$species_id) , ];
	res$vi$FN<-     res$vi$true[
			match( setdiff(res$vi$true$species_id, res$vi$pred$species_id), res$vi$true$species_id) , ];
	
	res$ba$true<-   true[true$superkingdom=="Bacteria",c("species","species_id","readn")];
	res$ba$pred<-   pred.ba[pred.ba$superkingdom=="Bacteria",c("species","species_id","readn")];
	res$ba$TP<-     res$ba$pred[
			match( intersect(res$ba$pred$species_id,res$ba$true$species_id), res$ba$pred$species_id) , ];
	res$ba$FP<-     res$ba$pred[
			match( setdiff(res$ba$pred$species_id,res$ba$true$species_id), res$ba$pred$species_id) , ];
	res$ba$FN<-     res$ba$true[
			match( setdiff(res$ba$true$species_id,res$ba$pred$species_id), res$ba$true$species_id) , ];
	
	table<- matrix(0,nrow=3*2,ncol=11);
	table<- data.frame(table);
	colnames(table)<- c("superkingdom","rank","pred","true","TP","FP","FN","Pr","Rc","F1","tail")
	king_list   <- c("Viruses","Bacteria");
	rank_list   <- c("species","genus","family");
	rankid_list <- c("species_id","genus_id","family_id");
	rowi        <- 1;
	for(i in 1:length(king_list)){
		for(j in 1:length(rank_list)){
			table[rowi,"superkingdom"]  <- king_list[i];
			table[rowi,"rank"]          <- rank_list[j];
			
			true_tmp<- unique( true[true$superkingdom==king_list[i],rankid_list[j] ] );
			pred_tmp = 0;
			if(king_list[i] == "Viruses"){
				pred_tmp<- unique( pred.vi[pred.vi$superkingdom==king_list[i],rankid_list[j] ] );
			}
			if(king_list[i] == "Bacteria"){
				pred_tmp<- unique( pred.ba[pred.ba$superkingdom==king_list[i],rankid_list[j] ] );
			}
			table[rowi,"pred"]  <- length(pred_tmp);
			table[rowi,"true"]  <- length(true_tmp);
			table[rowi,"TP"]    <- length(intersect(pred_tmp,true_tmp));
			table[rowi,"FP"]    <- length( setdiff(pred_tmp,true_tmp));
			table[rowi,"FN"]    <- length( setdiff(true_tmp,pred_tmp));
			if(king_list[i] == "Viruses"){
				table[rowi,"tail"]  <- tail_vi;
			}
			else if(king_list[i] == "Bacteria"){
				table[rowi,"tail"] <- tail_ba;   
			}
			rowi = rowi+1;
		}
	}
	# add "TOTAL" rows to summary table
	total<- table[1:3,];
	total[1,1]<- "TOTAL";
	total[1:3,3:10]<- table[1:3,3:10]+table[4:6,3:10];
	total[,11]<- NA
	table<- rbind(table,total);
	table[,"Pr"]<- table[,"TP"]/(table[,"TP"]+table[,"FP"]);
	table[,"Rc"]<- table[,"TP"]/(table[,"TP"]+table[,"FN"]);
	table$F1    <- (2*table$Pr*table$Rc)/(table$Pr+table$Rc);
	res$stats<- table;
	
	# if excel_file is given, write results to excel
	if(excel_file=="none"){ return(res);}
	
	wb<- createWorkbook();
	
	addWorksheet(wb, sheetName = "SUMMARY");
	addWorksheet(wb, sheetName = "vi_true");
	addWorksheet(wb, sheetName = "vi_pred");
	addWorksheet(wb, sheetName = "vi_TP");
	addWorksheet(wb, sheetName = "vi_FN");
	addWorksheet(wb, sheetName = "vi_FP");
	addWorksheet(wb, sheetName = "ba_true");
	addWorksheet(wb, sheetName = "ba_pred");
	addWorksheet(wb, sheetName = "ba_TP");
	addWorksheet(wb, sheetName = "ba_FN");
	addWorksheet(wb, sheetName = "ba_FP");
	
	writeDataTable(wb, sheet = "SUMMARY", x = res$stats, tableStyle = "None", withFilter = F);
	writeDataTable(wb, sheet = "vi_true", x = res$vi$true, tableStyle = "None", withFilter = F);
	writeDataTable(wb, sheet = "vi_pred", x = res$vi$pred, tableStyle = "None", withFilter = F);
	writeDataTable(wb, sheet = "vi_TP",   x = res$vi$TP, tableStyle = "None", withFilter = F);
	writeDataTable(wb, sheet = "vi_FN",   x = res$vi$FN, tableStyle = "None", withFilter = F);
	writeDataTable(wb, sheet = "vi_FP",   x = res$vi$FP, tableStyle = "None", withFilter = F);
	
	writeDataTable(wb, sheet = "ba_true", x = res$ba$true, tableStyle = "None", withFilter = F);
	writeDataTable(wb, sheet = "ba_pred", x = res$ba$pred, tableStyle = "None", withFilter = F);
	writeDataTable(wb, sheet = "ba_TP",   x = res$ba$TP, tableStyle = "None", withFilter = F);
	writeDataTable(wb, sheet = "ba_FN",   x = res$ba$FN, tableStyle = "None", withFilter = F);
	writeDataTable(wb, sheet = "ba_FP",   x = res$ba$FP, tableStyle = "None", withFilter = F);
	
	# set SUMMARY STYLE
	sty_pc  <- createStyle(numFmt = "0.0%");
	addStyle(wb, sheet = "SUMMARY", style = sty_pc,  rows = 1:(nrow(res$stats)+1), cols = c(8,9,10),gridExpand = T );
	
	# SET STYLE FOR THE REST
	sty_int <- createStyle(numFmt = "0");
	name_list <- names(wb);
	name_list <- setdiff(name_list,"SUMMARY");
	nrow_max  <- max(nrow(res$vi$true),nrow(res$vi$pred),nrow(res$ba$true),nrow(res$ba$pred));
	for(name in name_list){
		addStyle(wb, sheet = name, style = sty_int,  rows = 1:nrow_max, cols = c(2,3),gridExpand = T );
		setColWidths(wb, sheet = name, cols = 1:3, widths = "auto");
	}
	
	# Highlight false negatives/positives in vi_true/vi_pred spreadsheets
	sty_false<- createStyle(fgFill="#ff5050");
	rows.FN = match( setdiff(res$vi$true$species,res$vi$pred$species), res$vi$true$species) + 1;
	rows.FP = match( setdiff(res$vi$pred$species, res$vi$true$species), res$vi$pred$species) + 1;
	addStyle(wb, sheet="vi_true", sty_false, rows = rows.FN, cols= c(1:3), gridExpand=T);
	addStyle(wb, sheet="vi_pred", sty_false, rows = rows.FP, cols= c(1:3), gridExpand=T);
	
	# Highlight false negatives/positives in ba_true/ba_pred spreadsheets
	rows.FN = match( setdiff(res$ba$true$species,res$ba$pred$species), res$ba$true$species) + 1;
	rows.FP = match( setdiff(res$ba$pred$species, res$ba$true$species), res$ba$pred$species) + 1;
	addStyle(wb, sheet="ba_true", sty_false, rows = rows.FN, cols= c(1:3), gridExpand=T);
	addStyle(wb, sheet="ba_pred", sty_false, rows = rows.FP, cols= c(1:3), gridExpand=T);
	
	saveWorkbook(wb,excel_file,overwrite = T);
	return(res);
}


get_abund_tables<- function(abund_table, excel_file="none", 
		taxranks=c("species","genus","family","division"),
		taxgroups=c(), tail=0){
	#' @description  Reads raw abundance table: with readn [+contign] for taxa below species rank.
	#' And converts this to a abundance tables for Viruses/Bacteria at species/genus/family/division level
	#'
	#' @param abund_table file. Tab delimited file containing taxa abundancies. Format: taxid\t readn\t [contign]
	#' @param excel_file file. Results are printed here.
	#' @param taxranks vector. List of taxranks to report.
	#' @param taxgroups vector. List of taxgroups to report. These will be matched to the topmost taxrank and reported on separate sheets.
	#'							Setting taxgroups to c("all") will report all topmost taxranks.
	#' @param tail num. Number in percents [0,100], taxa with least abundance summing to this amount will be ignored
	#'  
	
	# Libraries + parameters
	library(reshape, quietly = T);
	library(openxlsx, quietly = T);
	toptaxrank		= taxranks[length(taxranks)];
	
	# START WORKING
	pred		= read.table(file=abund_table,sep="\t",header=T,comment.char="",stringsAsFactors = F, quote="", skipNul=T);
	for(rank in taxranks){
		if(!(rank %in% colnames(pred)) ){
			stop( sprintf("column %s not found in %s",rank,abund_table) );
		}
	}
	
	# convert empty/NA taxon labels to "unknown"
	for(rank in taxranks){
		pred[grep("^$",pred[[rank]]), rank] 	= "unknown";
		pred[is.na(pred[[rank]]), rank]		= "unknown";
	}
	# when reporting all available taxgroups
	if(length(taxgroups)== 0 ||
		(length(taxgroups)== 1 && taxgroups[1] =="all")){
		taxgroups = sort(unique( pred[[toptaxrank]]) );
	}
	
	# ADD PERCENTAGE COL
	readn_sum		= sum(pred$readn);
	pred$readn_pc	= pred$readn/readn_sum;
	
	# SUM OVER TAXRANKS
	res					= list();
	tmp					= NA; 
	if("contign" %in% colnames(pred)){
		tmp				= melt(pred[, c("readn","readn_pc","contign",taxranks) ], id=taxranks);
	}
	else{
		tmp				= melt(pred[, c("readn","readn_pc",taxranks) ], id=taxranks);
	}
	for(i in 1:length(taxranks)){
		rank			= taxranks[i];
		formula_tmp		= sprintf("%s~variable",paste(taxranks[i:length(taxranks)],sep = "+",collapse = "+"));
		res[[rank]]		= cast(tmp, formula_tmp, sum);
	}
	
	# ADD BPHAGE FIELD
	if("bphage" %in% colnames(pred)){
		for(rank in taxranks){
			res[[rank]]$bphage	= NA;
			for(i in 1:nrow(res[[rank]])){
				bphage_vals				= unique( pred$bphage[ pred[[rank]]==res[[rank]][i,rank] ]);
				res[[rank]]$bphage[i]	= paste(bphage_vals,sep="/",collapse = "/")
			}
		}
	}
	
	# RE-COUNT CONTIG NUMBER FOR EACH taxon
	#if("contign" %in% colnames(pred)){
	#	for(rank in taxranks){
	#		res[[rank]]$contign			= NA;
	#		for(i in 1:nrow(res[[rank]])){
	#			res[[rank]]$contign[i]	= sum( pred$contign[ pred[[rank]]==res[[rank]][i,rank] ] );
	#		}
	#	}
	#}
	
	# SORT
	for(rank in taxranks){
		res[[rank]]		= res[[rank]][ order(res[[rank]]$readn, decreasing=T),];
	}
	
	# CALC CUMSUM + DIVIDE CUMSUM INTO INTO CLASSES BY RANGE: 1:[0-95[, 2:[95-99[, 3:[99-100]
	cumsum_bounds		= c(1.0,0.99,0.95);
	cumsum_qlabels		= c(3,2,1);
	for(rank in taxranks){
		res[[rank]]$csum	=  (cumsum(res[[rank]]$readn) - res[[rank]]$readn)/ readn_sum;
		res[[rank]]$csumq	= cumsum_qlabels[1];
		for(i in 2:length(cumsum_bounds)){
			res[[rank]]$csumq[res[[rank]]$csum < cumsum_bounds[i]] = cumsum_qlabels[i];
		}
	}
	
	# REORDER COLS
	cols_pref		= c("readn","readn_pc","csum","csumq");
	if("contign" %in% colnames(pred)){
		cols_pref	= c(cols_pref,"contign");
	}
	if("bphage" %in% colnames(pred)){
		cols_pref	= c(cols_pref,"bphage");
	}
	for(i in 1:length(taxranks)){
		rank				= taxranks[i];
		taxranks_selected	= taxranks[i:length(taxranks)];
		res[[rank]]			= res[[rank]][,c(cols_pref,taxranks_selected)];
	}
	
	# RETRIEVE SUBGROUPS
	subgroup_list	= c();
	for(group in taxgroups){
		for(i in 1:(length(taxranks)-1)){
			rank			= taxranks[i];
			taxranks_sel	= taxranks[i:length(taxranks)];
			label			= sprintf("%s_%s",group,substr(rank,1,3));
			tmp				= NA;
			if("bphage" %in% colnames(res[[rank]]) & (group=="Viruses")){
				tmp			= res[[rank]][ res[[rank]][[toptaxrank]] == group,c(cols_pref,taxranks_sel) ];
			}
			else{
				tmp			= res[[rank]][ res[[rank]][[toptaxrank]] == group,c(setdiff(cols_pref,"bphage"),taxranks_sel) ];
			}
			if(nrow(tmp)>0){
				res[[label]]= tmp;
				subgroup_list	= c(subgroup_list,label);
			}
		}
	}
	
	# ADD ROW NUMS FOR CLARITY (comment out to preserv row numbers from the original pred_file)
	for(name in subgroup_list){
		if(nrow(res[[name]]) > 0){
			rownames(res[[name]])<- c(1:nrow(res[[name]]));
		}
	}
	
	# WRITE TO EXCEL FILE  (IF GIVEN)
	if(excel_file=="none"){ return(res);}
	wb			= createWorkbook();
	sty_int <- createStyle(numFmt = "0");
	sty_pc  <- createStyle(numFmt = "0.0%");
	sty_tail <- createStyle(fontColour = "#808080");
	
	for(name in subgroup_list){
		if(nrow(res[[name]]) > 0){
			addWorksheet(wb,sheetName=name);
			writeDataTable(wb,sheet=name,x= res[[name]],tableStyle = "TableStyleMedium4");
			# STYLE:
			M<- nrow(res[[name]]) + 1;
			N<- ncol(res[[name]]);
			tail_rows<- which(res[[name]]$csum >= ((100-tail)/100.0) ) + 1;
			addStyle(wb, sheet = name, style = sty_int, rows = 2:M, cols = 1);
			addStyle(wb, sheet = name, style = sty_pc,  rows = 2:M, cols = 2);
			addStyle(wb, sheet = name, style = sty_pc,  rows = 2:M, cols = 3);
			setColWidths(wb, sheet = name, cols = 1:N, widths = "auto");
			
			if( (name!= "sking") && (length(tail_rows) > 0)){ # cutting tail does not work for taxa at highest rank
				addStyle(wb, sheet = name, style = sty_tail,rows = tail_rows, cols = c(1:N),gridExpand = T,stack = T);
			}
		}
	}
	saveWorkbook(wb,excel_file,overwrite = T);
	return(res);
}

# 
# 
# Parameters:
# annot_file       : path to the annotation file
# excel_file       : Write generated annotation tables as data sheets to this file (overwrite if exists).
# toptaxrank	   : Top taxonomy rank for subdividing results into separate sheets
#
#' @description  GENERATES SEQUENCE ANNOTATION TABLES
#'
#' @param annot_file file. Tab delimited file containing taxa abundancies. Format: taxid\t readn\t [contign]
#' @param excel_file file. Write generated annotation tables as data sheets to this file (overwrite if exists).
#' @param toptaxrank str.  Taxrank used to split results into tables/datasheers.
#' @param taxgroups vector. List of taxgroups to report. These will be matched to the topmost taxrank and reported on separate sheets.
#'							Setting taxgroups to c("all") will report all topmost taxranks.
#'  
get_annot_tables<- function(annot_file,excel_file="none",
		toptaxrank="division",
		taxgroups=c() ){
	# Libraries and settings
	library(openxlsx, quietly = T);
	
	
	# START WORKING
	pred<- read.table(file=annot_file,sep="\t",header=T,comment.char="",stringsAsFactors = F,quote="", skipNul=T);
	if(!(toptaxrank %in% colnames(pred)) ){
		stop( sprintf("column %s not found in %s",toptaxrank,annot_file) );
	}
	
	# convert empty/NA taxon labels to "unknown"
	rank= toptaxrank;
	pred[grep("^$",pred[[rank]]), rank] 	= "unknown";
	pred[is.na(pred[[rank]]), rank]			= "unknown";
	
	# when reporting all available taxgroups
	if(length(taxgroups)==0 ||
		(length(taxgroups)==1 && taxgroups[1]=="all") ){
		taxgroups = sort(unique( pred[[toptaxrank]]) );
	}
	
	# delete *_id cols
	ind		= grep("_id",colnames(pred));
	if(length(ind)>0){
		pred	= pred[,-ind];
	}
	# select phage-col
	phage_colname	= NA;
	tmp				= grep("phage",colnames(pred),ignore.case = T,value=T);
	if(length(tmp)>0){
		phage_colname= tmp[1];
	}
	
	# SELECT SUBSETS: Viruses/Bacteria/..
	qseqid_colname = "qseqid";
	if("contig" %in% colnames(pred)){
		qseqid_colname = "contig";
	}
	
	res<- list();
	for(group in taxgroups){
		tmp = NA;
		if(group == "Viruses"){
			if(!is.na(phage_colname)){
				tmp = pred[ pred[[toptaxrank]] == group & pred[[phage_colname]]!="yes", ];
			}
			else{
				tmp = pred[ pred[[toptaxrank]] == group, ];
			}
		}
		else if(group == "Phages"){
			if(!is.na(phage_colname)){
				tmp = pred[ !is.na(pred[[phage_colname]]) & pred[[phage_colname]]=="yes", ];
			}
			else{
				tmp = pred[ pred[[toptaxrank]] == group, ];
			}
		}
		else{
			tmp = pred[ pred[[toptaxrank]] == group, setdiff(colnames(pred),c(phage_colname,"host.source","genome.composition"))];
		}
		if(nrow(tmp)>0){
			res[[group]]	= tmp;
		}
	}
	
	# SORT
	qseqlen_colname    = "qseqlen";
	if("clen" %in% colnames(pred)){
		qseqlen_colname = "clen";
	}
	for(group in names(res)){
		res[[group]]	= res[[group]][ order(res[[group]][,qseqlen_colname], res[[group]][,qseqid_colname], res[[group]]$bitscore, decreasing=T), ];
	}
	
	# PRINT
	if(excel_file=="none"){ return(res);}
	
	wb<- createWorkbook();
	for(name in names(res)){
		addWorksheet(wb,sheetName = name);
		writeDataTable(wb,sheet = name,x = res[[name]], tableStyle = "TableStyleMedium4");
		setColWidths(wb, sheet = name, cols = 1:ncol(res[[name]]), widths = "auto");
	}
	saveWorkbook(wb,excel_file,overwrite = T);
	return(res);
}

# 
# PROJECT:	LAZYPIPE
# SUBJECT:	FUNCTIONS FOR LAZYAPP REPORTING
# DATE: 	2024/02/26


read_hists<- function(read1, read2=NA, outdir, outpref="read.hist",
		title="Read length", ylab="Number of reads", xlab="Length (bp)",
		width=5.2, res=300, pointsize=9, binwidth=10){
	if(length(grep(".(len|length)$",read1,perl=T,ignore.case=T)) == 0){
		stop("invalid input\n", call.=FALSE);
	}
	library(ggplot2);
	library(cowplot);
	data1	= read.csv(file=read1, sep="", header=FALSE);
	data2	= NA;
	colnames(data1) = c("length");
	if(!is.na(read2)){
		data2= read.csv(file=read2, sep="", header=FALSE);
		colnames(data2)= c("length");
	}
	col.grey	= "#e9ecef";
	col.blue	= "#1e81b0";
	col.dblue	= "#063970";
	title1	= "";
	title2	= "";
	p1		= ggplot(data1, aes(x=length)) + 
			geom_histogram(binwidth = binwidth, fill=col.blue, color=col.dblue) +
			labs(title=title1) +
			xlab(label = xlab)+
			ylab(label = ylab)+
			theme(plot.title = element_text(size=12, hjust = 0.5));
	p2	= NA;
	if(!is.na(read2)){
		title1 = "Forward reads";
		title2 = "Reverse reads";
		p2		= ggplot(data2, aes(x=length)) + 
				geom_histogram(binwidth = binwidth, fill=col.blue, color=col.dblue) +
				labs(title=title2) +
				xlab(label = xlab)+
				ylab(label = ylab)+
				theme(plot.title = element_text(size=12, hjust = 0.5));
	}
	# save histogram(s)
	filename	= sprintf("%s.png",outpref);
	if(is.na(read2)){
		ggsave(filename,p1,png,outdir,width=width*res,height=width*res,units ="px",dpi =res);
	}
	else{
		panel = plot_grid(p1, p2, labels = c("A", "B"), nrow = 1, ncol=2);
		ggsave(filename,panel,png,outdir,width=width*res*2,height=width*res,units ="px",dpi =res);
	}
}
#
# Plots a number of histograms visualizing contig length data
#
contig_hists<- function(file_in, outdir, outpref="contig.hist",
		title="Contig length", ylab="Number of contigs", xlab="Length (bp)",
		width=5.2, res=300, pointsize=9, binwidth=10){
	if(length(grep(".(len|length)$",file_in,perl=T,ignore.case=T)) == 0){
		stop("invalid input\n", call.=FALSE);
	}
	library(ggplot2);
	library(cowplot);
	data	= read.csv(file=file_in, sep="", header=FALSE);
	colnames(data) = c("length");
	# break data by quantines
	len.min	= min(data$length);
	q1		= round(quantile(data$length, probs=c(0.01)));
	q5 		= round(quantile(data$length, probs=c(0.05)))
	q25		= round(quantile(data$length, probs=c(0.25)))
	q50 	= round(quantile(data$length, probs=c(0.50)))
	q75 	= round(quantile(data$length, probs=c(0.75)))
	q95		= round(quantile(data$length, probs=c(0.95)))
	q99 	= round(quantile(data$length, probs=c(0.99)))
	len.max	= max(data$length);
	grey	= "#e9ecef";
	blue	= "#1e81b0";
	dblue	= "#063970";
	mdata	= data.frame(
			from =c(len.min,len.min,len.min,q1,q5,q95,q99),
			to =c(len.max,q1,q5,q99,q95,len.max,len.max),
			#bins=c(100,10,10,50,50,10,10),
			name= c('q0-100','q0-1','q0-5','q1-99','q5-95','q95-100','q99-100'),
			range=c('','Q0..Q1', 'Q0..Q5','Q1..99','Q5..Q95','Q95..Q100','Q99..Q100'),
			fill= rep(blue,7),
			color= c(dblue,dblue,dblue,dblue,dblue,dblue,dblue));
	
	plot_list	= list();
	for(i in 1:nrow(mdata)){
		data.tmp = data[data$length>= mdata$from[i] & data$length<=mdata$to[i],1,drop=F];
		title	= sprintf("%s (%s to %s bp)",
				mdata$range[i],
				format(mdata$from[i],big.mark=","),
				format(mdata$to[i],big.mark=","));
		plot_list[[ mdata$name[i] ]] = 
				ggplot(data.tmp, aes(x=length)) + 
				geom_histogram(binwidth = binwidth, fill=mdata$fill[i], color=mdata$color[i],) +
				labs(title=title) +
				xlab(label = xlab)+
				ylab(label = ylab)+
				theme(plot.title = element_text(size=12, hjust = 0.5));
		#ggsave(filename,p,png,outdir,width=width*res,height=width*res,units ="px",dpi =res);
	}
	
	# save histograms combining some plots into multipanels
	# q0-100 hist
	filename	= sprintf("%s.png",outpref);
	ggsave(filename,plot_list[['q0-100']],png,outdir,width=width*res*2,height=width*res,units ="px",dpi =res);
	
	# (A) q0-1 (B) q1-99 (C) q99-100 panel plot
	filename	= sprintf("%s.%s.png",outpref,'q1-99-100');
	panel		= plot_grid(plot_list[['q0-1']], plot_list[['q1-99']], plot_list[['q99-100']],
			labels = c("A", "B", "C"), nrow = 1, ncol=3);
	ggsave(filename,panel,png,outdir,width=width*res*3,height=width*res,units ="px",dpi =res);
	# (A) q0-5 (B) q5-95 (C) q95-100 panel plot
	filename	= sprintf("%s.%s.png",outpref,'q5-95-100');
	panel		= plot_grid(plot_list[['q0-5']], plot_list[['q5-95']], plot_list[['q95-100']],
			labels = c("A", "B", "C"), nrow = 1, ncol=3);
	ggsave(filename,panel,png,outdir,width=width*res*3,height=width*res,units ="px",dpi =res);	
}
