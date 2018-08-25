#!/usr/bin/env python

import argparse
from collections import defaultdict
import math
import pandas as pd
import re
import random 
from shutil import copyfile
import subprocess
import sys
import os.path

##########################################################################################################
##########################################################################################################
## Author: Daria Gavriouchkina (c)                                                                      ##
##                                                                                                      ##
## Given a subread featureCounts generated count table for multiple samples with replicates, calculates ##
## FPKM/RPKM values and selects samples expressing a gene for which Ensembl geneID is provided based on ##
## a threshold value. Samples are split into 2 categories - those expressing and those not expressing   ##
## the gene of interest. Randomly a given number of cells from pool of single cell data and determine   ##
## how reproducible pooling is.                                                                         ##
##########################################################################################################
##########################################################################################################


def make_list_from_column_file(filename):
	glist=[]
	for line in open(filename):
		glist.append(line.rstrip())
	print len(glist)

	return glist

def pick_randomly(list_of_cells, hm_cells, seed):
 	chosen=[]
 	random.seed(seed)
 	for i in range(0, int(hm_cells)):
 		cells2pick=random.choice(list_of_cells)
 		chosen.append(cells2pick)
 		try:
	 		cells_left = list_of_cells
	 		cells_left.remove(cells2pick)
 		except ValueError:
 			pass
 	return chosen, cells_left


def picking_cells(list_of_cells2, hm_cells2pick, nb_repl, seed):
	selected=[]
	left=[]
	while nb_repl > 0: 
		#print '\treplicate no=', nb_repl
		#print '\tTotal # of cells from which to pick=', len(list_of_cells2),'\n'
		(chosen_cells, cells_left)=pick_randomly(list_of_cells2, hm_cells2pick, seed)
		selected.append(chosen_cells)
		list_of_cells2=cells_left
		nb_repl -= 1 
	list_of_cells2 =None
	cells_left =None
	return selected

def write_list_to_file(gps_of_cells, filename):
	out=open(filename, 'w')
	for gpcell in gps_of_cells: 
		out.write(','.join(gpcell)+'\n')

def fpkm_calc(filename, geneid, limit):
	fpkm_r = open('fpkm.R', 'w')
	fpkm_r.write('''countToFpkm <- function(counts, effLen)
{
	N <- sum(counts)
	exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}
''')
	fpkm_r.write('cts <- read.table("'+ filename+'", sep="\\t", row.names=1, h=T)\n')
	fpkm_r.write('cts$Chr <- NULL\ncts$Start <- NULL\ncts$End <- NULL\ncts$Strand <- NULL\n')
	fpkm_r.write('old_ncol <- ncol(cts)\n')

	cmd2="sed -n '2p'  {0} > {1}".format(filename,"all_cell_names.txt")
	os.system(cmd2)

	with open('all_cell_names.txt', 'r') as infile:
		colnames = infile.read().rstrip().split('\t') 

	removed_cnames=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']
	for cname in colnames: 
		if cname not in removed_cnames:
			fpkm_r.write('cts$'+cname+'_fpkm <- countToFpkm(cts$'+cname+', cts$Length)\n') 
	fpkm_r.write('start <- old_ncol +1\n')
	fpkm_r.write('cts.fpkm <- cts[,c(start:ncol(cts))]\n')
	fpkm_r.write('write.table(as.data.frame(cts.fpkm), file="'+'FPKM_'+filename+ '", sep="\\t", quote=F)\n') # cts.fpkm colnames=cell names, rownames=genes
	fpkm_r.write('tcts.fpkm <- as.data.frame(t(cts.fpkm))\n')   # transpose so colnames=genes, rownames=cells
	fpkm_r.write('tcts.fpkm.pos_list <- rownames(tcts.fpkm[tcts.fpkm$'+geneid+' > '+ str(limit)+',])\n')
	fpkm_r.write('tcts.fpkm.neg_list <- rownames(tcts.fpkm[!rownames(tcts.fpkm ) %in% tcts.fpkm.pos_list,]) \n')
	fpkm_r.write('write.table(as.data.frame(tcts.fpkm.pos_list), file="'+geneid+'_expressing_cell_ids.txt", quote=F, row.names=F, col.names=F)\n')
	fpkm_r.write('write.table(as.data.frame(tcts.fpkm.neg_list), file="'+geneid+'_NOT_expressing_cell_ids.txt", quote=F, row.names=F, col.names=F)\n')
	
	return colnames

def merge_and_run_DESeq2(fname, nb_replics, nb_cells, iteration_nb):
	print iteration_nb
	merge_r =open('mergeCells'+'_'+str(iteration_nb)+'.R', 'w')
	merge_r.write('library(DESeq2)\nlibrary(fdrtool)\n')
	merge_r.write('cts <- read.table("'+ fname+'", sep="\\t", row.names=1, h=T)\n')
	merge_r.write('cts$Chr <- NULL\n')
	merge_r.write('cts$Start <- NULL\n')
	merge_r.write('cts$End <- NULL\n')
	merge_r.write('cts$Strand <- NULL\n')
	merge_r.write('cts$Length <- NULL\n')

	merge_r.write('expr <- read.table("'+str(iteration_nb) + '_expr.txt", sep=",", h=F)\n')
	merge_r.write('texpr <- as.data.frame(t(expr))\n')
	
	merge_r.write('noexpr <- read.table("'+str(iteration_nb)+'_noexpr.txt", sep=",", h=F)\n')
	merge_r.write('tnoexpr <- as.data.frame(t(noexpr))\n')
	
	expr_listV=[]
	noexpr_listV=[]
	
	for i in range(1, nb_replics+1):
		
		expr_listV.append('cts$exprSum_repl_'+str(i)) 
		noexpr_listV.append('cts$noexprSum_repl_'+str(i)) 
		
		merge_r.write('cts$exprSum_repl_'+str(i)+' <- rowSums(cts[,colnames(cts) %in% texpr$V'+str(i)+'])\n')
		merge_r.write('cts$noexprSum_repl_'+str(i)+' <- rowSums(cts[,colnames(cts) %in% tnoexpr$V'+str(i)+'])\n')
	
	merge_r.write("cts2 <- data.frame(rownames(cts)," + ",".join(expr_listV)+ "," + ",".join(noexpr_listV) + ")\n")
	merge_r.write("rownames(cts2) <- cts2$rownames.cts.\n")
	merge_r.write("cts2$rownames.cts. <- NULL\n")
	
	# Diff Exp 
	
	merge_r.write('de.dsg <- data.frame(row.names=colnames(cts2), type=as.factor(c("positive"'+', "positive" '*(nb_replics-1) +  ', "negative"'*(nb_replics-1) +',"negative")))\n')
	merge_r.write('de.dsq <- DESeqDataSetFromMatrix(countData=cts2, colData=de.dsg, design=~type)\n')
	merge_r.write('de.dsq <- DESeq(de.dsq)\n')
	merge_r.write('old_de.res <- results(de.dsq, contrast=c("type", "positive",  "negative"))\n')
	merge_r.write('old_de.res_up <- old_de.res[old_de.res$log2FoldChange > 0 & old_de.res$padj < 0.05 & !is.na(old_de.res$padj), ]\n')
	merge_r.write('old_de.res_dn <- old_de.res[old_de.res$log2FoldChange < 0 & old_de.res$padj < 0.05 & !is.na(old_de.res$padj), ]\n')

	# Write to file results of diff exp
	
	merge_r.write('write.table(as.data.frame(old_de.res_up), file="up_'+str(iteration_nb)+'.txt", sep="\\t", quote=F)\n')
	merge_r.write('write.table(as.data.frame(old_de.res_dn), file="dn_'+str(iteration_nb)+'.txt", sep="\\t", quote=F)\n')

	merge_r.write('pdf("QC_Histogram_PVal_Padj_'+str(nb_cells)+'_Cells_Iteration'+str(iteration_nb)+'.pdf")\n')
	merge_r.write('hist(old_de.res$pvalue, col = "dodgerblue", main ="Histogram PVal'+ str(nb_cells)+'\\nCells Iteration' + str(iteration_nb)+'", xlab = "p-values")\n')
	merge_r.write('hist(old_de.res$padj, col = "dodgerblue", main ="Histogram PAdj'+ str(nb_cells)+'\\nCells Iteration ' + str(iteration_nb)+'", xlab = "Adjusted p-values")\n')
	merge_r.write('dev.off()\n')
	
	merge_r.write('new_de.res <- old_de.res\n')
	merge_r.write('new_de.res <- new_de.res[ !is.na(new_de.res$padj), ]\n')
	merge_r.write('new_de.res <- new_de.res[ !is.na(new_de.res$pvalue), ]\n')
	merge_r.write('new_de.res <- new_de.res[, -which(names(new_de.res) == "padj")]\n')
	merge_r.write('FDR.test.res <- fdrtool(new_de.res$stat, statistic= "normal", plot = T)\n')

	merge_r.write('if(FDR.test.res$param[1, "sd"] != 1)\n')
	merge_r.write('\tnew_de.res[,"padj"]  <- p.adjust(FDR.test.res$pval, method = "BH")\n')
	merge_r.write('\tpdf("Corrected_QC_Histogram_PVal_Padj_'+str(nb_cells)+'_Cells_Iteration'+str(iteration_nb)+'.pdf")\n')
	merge_r.write('\tpar(mar=c(4,4,4,4))\n')
	merge_r.write('\thist(new_de.res$pval, col = "dodgerblue", main = "Histogram PVal'+str(nb_cells)+'\\nCells Iteration'+str(iteration_nb)+'",  xlab = "CORRECTED p-values")\n')
	merge_r.write('\tdev.off()\n')
	merge_r.write('\tnew_de.res_up <- new_de.res[new_de.res$log2FoldChange > 0 & new_de.res$padj < 0.05 & !is.na(new_de.res$padj), ]\n')
	merge_r.write('\tnew_de.res_dn <- new_de.res[new_de.res$log2FoldChange < 0 & new_de.res$padj < 0.05 & !is.na(new_de.res$padj), ]\n')
	
	merge_r.write('\twrite.table(as.data.frame(new_de.res_up), file="PVal_Corr_up.txt", sep="\\t", quote=F, col.names=F)\n')
	merge_r.write('\twrite.table(as.data.frame(new_de.res_dn), file="PVal_Corr_dn.txt", sep="\\t", quote=F, col.names=F)\n')
	
	
	
if __name__ == "__main__":	
	
	# Parse Arguments 
	parser = argparse.ArgumentParser(description="Attempt to parse single cell count table generated by subread FeatureCounts based on FPKM values and create multiple cell pools with given number of iterations")
	
	parser.add_argument('-f', '--filename', type=str, dest='filename', action='store', required=True, help='filename for single cell count table produced by subread FeatureCounts')
	parser.add_argument('-gi', '--gene_id', type=str, dest='geneid',  action='store', required=True, help='Ensembl Geneid based on which to parse count table')
	parser.add_argument('-o', '--output_dir', type=str, dest='outdir',  action='store', required=True, help='name of directory in which to store analyses')
	parser.add_argument('-c', '--no_cells',type=int, dest='no_cells',  action='store', help='number of cells to pool (default=50)', default=50)
	parser.add_argument('-l', '--limit',  type=int, dest='limit', default=1,  action='store', help="FPKM cutoff to consider gene expressed (default is 1 FPKMs)")
	parser.add_argument('-i', '--no_iterations', type=int, dest='iterations', default=100, action='store', help='Number of times to resample cells and run DESeq2')
	parser.add_argument('-v', '--verbose', action="store_true", help='verbose output')
	parser.add_argument('-V', '--version', action='version', version='%(prog)s (version 0.1)')	

	args = parser.parse_args()
	#print args
	# eg.
	# Namespace(filename='STAR_mapped_zv10_single_cell_satija_counts', geneid='ENSDARG00000021032', limit=2, no_cells=50, verbose=True)
	print '\n\n\n'
	print '#######################################################################################'
	print '####                                                                                ###'
	print '####                       Single_cell_iterative_pooling.py                         ###'
	print '####                                                                                ###'
	print '#######################################################################################\n\n\n'

	print '#######################################################################################\n'

	###########################################
	##########  Check arguments ... ###########
	###########################################
	
	if args.verbose ==True:
		print '\tLooking for file "', args.filename, '" ....'
	if os.path.isfile(args.filename): 
		if args.verbose ==True:
			print '\t\tfile ',  args.filename, ' exists ! \n'
	else: 
		print '\tno file called ', args.filename, ' exists in ', os.getcwd()
		sys.exit("\tTry again!")
	
	if args.verbose ==True:
		print '\tLooking for gene "', args.geneid, ' in file', args.filename, '" ...'
	cmd = "tail -n +2 {0} | cut -f1  >  {1}".format(args.filename, 'all_gene_names.txt')
	os.system(cmd)
	glist=make_list_from_column_file('all_gene_names.txt')
	
	if args.geneid in glist: 
		if args.verbose ==True:
			print '\tGene "', args.geneid, '" successfully found in file "', args.filename, '" !\n'
	else: 
		print '\tGene ', args.geneid, ' was NOT found in file', args.filename
		sys.exit("Try again!")
	
	directory= '/'.join([os.getcwd(), args.outdir])
	if os.path.isdir(directory):
		print 'The output directory',  args.outdir, ' already exists. Please come up with a new name.'
		sys.exit("\tTry again!")
	if (bool(re.compile('^[A-Za-z0-9_\.]+$').match(args.outdir.rstrip())) ==False) :
		print 'Strange character in the directory name you provided: ', args.outdir
		sys.exit("Try again!")
		
	if args.iterations <= 1: 
		print '\tNumber of iterations must be more than 1! '
		sys.exit("\tTry again!")
	
	if args.verbose ==True:	
		print '\tOther parameters used: '
		print '\t########################'

	if args.limit < 0: 
		print '\tCutoff for expression must be a positive number'
		sys.exit("\tTry again!")
		
	if args.verbose ==True:
		print '\tCutoff for expression = ', args.limit
	
	if args.no_cells <= 1 : 
		print '\tNumber of cells to pool in each group must be > 1!'
		sys.exit("Try again!")
	
	if args.verbose ==True:
		print '\tNumber of cells to pool together', args.no_cells	
	
	if args.verbose ==True:
		print '\tVerbose mode on!\n'
	
	
	##################################################################################################
	### Creating output directory and moving count table file into it                             ####
	##################################################################################################
	
	
	if not os.path.exists(directory):
		os.makedirs(directory)
	
	src='/'.join([os.getcwd(), args.filename])
	#print src
	dest='/'.join([directory, args.filename])
	#print dest
	copyfile(src, dest)
	os.chdir(directory)
	
	
	
	print '\n\n\n'
	print '\t#############################'
	print '\t# Calculatig FPKM values  : # '
	print '\t#############################\n'

	
	# Calculate FPKMs based on count table, write to file and determine which cells express gene of interest 
	
	colnames = fpkm_calc(args.filename, args.geneid, args.limit)
	
	cmd3="Rscript fpkm.R"
	os.system(cmd3)
	print '\tDone making FPKMs!\n '
	
	expr = make_list_from_column_file(args.geneid+'_expressing_cell_ids.txt')
	expr2 =  ['_'.join(a.split('_')[:-1]) for a in expr]
	ncol_expr=len(expr)
	print 'EXPR', len(expr2), type(expr2)
	
		
	noexpr= make_list_from_column_file(args.geneid+'_NOT_expressing_cell_ids.txt')
	noexpr2 =  ['_'.join(b.split('_')[:-1]) for b in noexpr]
	ncol_noexpr = len(noexpr)
	
	print 'NOT EXPR', len(noexpr2), type(expr2)

	if args.verbose ==True: 
	 	print '\tParsing FPKM table by ', args.geneid
	 	print '\t########################################\n'
	 	print '\tNo cells ARE expressing ', args.geneid, 'at ',  str(args.limit), '= ', ncol_expr
	 	print '\tNo cells NOT expressing ', args.geneid, 'at ',  str(args.limit), '= ', ncol_noexpr
	 	print '\t-----------------------------------------------------------------'
	 	print '\tTotal number of cells= ', len(colnames)-6, '\n'

	#Check that number of cells requested is possible with this gene 

	print '\tChecking that requested cells pools are possisble : '
	print '\t###################################################\n'
	
	print "ncol_expr/args.no_cells=", ncol_expr/args.no_cells, type(ncol_expr/args.no_cells)
	print ncol_expr/args.no_cells <= 2 
	
	if ncol_expr/args.no_cells < 2 or ncol_noexpr/args.no_cells < 2 :
		print '\tNot enough cells to split: ',str(ncol_expr),'/' ,str(len(colnames)-6), 'expressing ', args.geneid
	 	print '\tCannot split ',str(ncol_expr), 'into ',  args.no_cells
	 	sys.exit("\tTry again!")
	else: 
		print '\tOK - will do ... \n'

	no_repl= min(math.trunc(ncol_expr/args.no_cells), math.trunc(ncol_noexpr/args.no_cells))


	if args.verbose ==True: 
		print '\tDetermining number of replicates '
		print '\t##################################'
	
		print '\tRatio ARE expressing cells : requested cells = ', ncol_expr, ': ', args.no_cells, '=', math.trunc(ncol_expr/args.no_cells)
		print '\tRatio NOT expressing cells : requested cells =', ncol_noexpr, ': ', args.no_cells, '=', math.trunc(ncol_noexpr/args.no_cells)
	
	
		print '\tNumber of replicates ==>', no_repl
		print '\t--------------------------\n'
	
	u_gene_boot_dict=defaultdict(int)
	u_gene_pvallogFC_dict=defaultdict(list)

	d_gene_boot_dict=defaultdict(int)
	d_gene_pvallogFC_dict=defaultdict(list)

	ua_gene_boot_dict=defaultdict(int)
	ua_gene_pvallogFC_dict=defaultdict(list)

	da_gene_boot_dict=defaultdict(int)
	da_gene_pvallogFC_dict=defaultdict(list)

	summary=open(str(args.geneid)+'_over_'+str(args.limit)+'_FPKM_in_'+str(args.no_cells)+'_Cells_'+str(args.iterations)+'iterations.log', 'w')
	
	summary.write('\t'.join(['Iteration#','Seed','PositiveCellsPooled', 'NegativeCellsPooled', '#UpregulatedGenes', '#DwnRegualtedGenes', 'AfterCorr#UpregulatedGenes', 'AfterCorr#DwnRegulatedGenes' ])+'\n')
	summary.write('------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
	
	
	###########################
	## Start of iteration    ##
	###########################
	#print 'BEFORE ITERATION\texpr, noexpr2', len(expr2), len(noexpr2), type(expr2), type(noexpr2), id(expr2), id(noexpr2)

	list_of_seeds =[]

	for it in range (1,(int(args.iterations)+1)): 
		
		#print '\tIteration #', it
		rexpr2  = expr2[:]
		rnoexpr2  = noexpr2[:]
		#print 'AFTER START ITERATION\texpr, noexpr2\t', len(expr2), len(noexpr2), type(expr2), type(noexpr2), id(expr2), id(noexpr2)
		#print 'AFTER [:] COPY\tRexpr, Rnoexpr2\t', len(rexpr2), len(rnoexpr2), type(rexpr2), type(rnoexpr2), id(rexpr2), id(rnoexpr2)

		# randomly generate a seed 

		seed = random.randint(0, sys.maxint)
		list_of_seeds.append(seed)
		
		if args.verbose== True:
			print '\tIteration #', it
			print '\t###################\n'

		if args.verbose== True:
			print '\tseed is ', seed, '\n'
			print '\tPicking expressing cells'
			print '\t########################\n'
		
		selected_expr_cells = picking_cells(rexpr2, args.no_cells, no_repl, seed)
		
		str_sel_expr= ";".join(",".join(map(str,gp_sel_cells)) for gp_sel_cells in selected_expr_cells)
		
		#print 'AFTER PICKING EXPRESSED:\tRexpr, Rnoexpr2\t', len(rexpr2), len(rnoexpr2), type(rexpr2), type(rnoexpr2), id(rexpr2), id(rnoexpr2)


		if args.verbose== True:
			print '\n'
			print '\tPicking NON-expressing cells'
			print '\t########################\n'

		selected_noexpr_cells = picking_cells(rnoexpr2, args.no_cells, no_repl, seed)	
		str_sel_noexpr= ";".join(",".join(map(str,gpno_sel_cells)) for gpno_sel_cells in selected_noexpr_cells)


		#print 'AFTER PICKING NOt EXPRESSED\tRexpr, Rnoexpr2\t', len(rexpr2), len(rnoexpr2), type(rexpr2), type(rnoexpr2), id(rexpr2), id(rnoexpr2)
		
		write_list_to_file(selected_expr_cells, str(it)+'_expr.txt')
		write_list_to_file(selected_noexpr_cells, str(it)+'_noexpr.txt')
		
		#print 'AFTER WRITING SELECTED CELLS TO FILE\tRexpr, Rnoexpr2\t', len(rexpr2), len(rnoexpr2), type(rexpr2), type(rnoexpr2), id(rexpr2), id(rnoexpr2)

		
		if args.verbose ==True:
			print '\tMerging cells and running DESeq2...\n'
		
		merge_and_run_DESeq2(args.filename, no_repl, args.no_cells, it)
		
		#print 'AFTER MERGING AND RUNNING DESEQ2\tRexpr, Rnoexpr2\t', len(rexpr2), len(rnoexpr2), type(rexpr2), type(rnoexpr2), id(rexpr2), id(rnoexpr2)
		
		cmd4="Rscript mergeCells"+"_"+str(it)+".R"
		os.system(cmd4)
		
		if args.verbose ==True:
			print '\tFinished doing that ... \n'
			print '\tNow going through DESeq2 results : \n'
		
		uct=0
		for ul in open('up_'+str(it)+'.txt'): 
			if not ul.startswith('baseMean'):
				uct+=1
				u=ul.rstrip().split('\t')
				u_gene_boot_dict[u[0]]+=1
				u_gene_pvallogFC_dict[u[0]].append([str(seed),u[2], u[6]])

		
		dct=0
		for dl in open('dn_'+str(it)+'.txt'): 
			if not dl.startswith('baseMean'):
				dct+=1
				d=dl.rstrip().split('\t')
				d_gene_boot_dict[d[0]]+=1
				d_gene_pvallogFC_dict[d[0]].append([str(seed), d[2], d[6]])
		
		uact=0
		for ual in open('PVal_Corr_up.txt'):
			if not ual.startswith('baseMean'):	
				uact+=1
				ua=ual.rstrip().split('\t')
				ua_gene_boot_dict[ua[0]]+=1
				ua_gene_pvallogFC_dict[ua[0]].append([str(seed), ua[2],ua[6]])

		
		dact=0
		for dal in open('PVal_Corr_dn.txt'):
			if not dal.startswith('baseMean'):
				dact+=1
				da=dal.rstrip().split('\t')
				da_gene_boot_dict[da[0]]+=1
				da_gene_pvallogFC_dict[da[0]].append([str(seed), da[2],da[6]])
		
		
		#print 'AFTER MAKING DICTIONARIES\tRexpr, Rnoexpr2', len(rexpr2), len(rnoexpr2), type(rexpr2), type(rnoexpr2), id(rexpr2), id(rnoexpr2)

		if args.verbose ==True:
			print '\t\tUp\tDn\tAdjUp\tAdjDwn\n'
			print '------------------------------------------------'
			print '\tit1\t',uct,'\t',dct,'\t',uact,'\t', dact, '\n'
		
		summary.write('\t'.join([str(it), str(seed), str_sel_expr, str_sel_noexpr, str(uct), str(dct), str(uact), str(dact) ])+'\n')	
				
	#Only bootstraps
	
	up_out_boot= open(str(args.geneid)+'_over_'+str(args.limit)+'_FPKM_in_'+str(args.no_cells)+'_Cells_'+str(args.iterations)+'iterations_up_bt.txt', 'w')
	
	for ku_gene_boot,vu_gene_boot in u_gene_boot_dict.iteritems():
		up_out_boot.write('\t'.join([ku_gene_boot,str(vu_gene_boot)])+'\n')
		
	dn_out_boot=open(str(args.geneid)+'_over_'+str(args.limit)+'_FPKM_in_'+str(args.no_cells)+'_Cells_'+str(args.iterations)+'iterations_dn_bt.txt', 'w')
	
	for kd_gene_boot,vd_gene_boot in d_gene_boot_dict.iteritems():
		dn_out_boot.write('\t'.join([kd_gene_boot,str(vd_gene_boot)])+'\n')


	up_adj_boot=open(str(args.geneid)+'_over_'+str(args.limit)+'_FPKM_in_'+str(args.no_cells)+'_Cells_'+str(args.iterations)+'iterations_up_adj_bt.txt', 'w')
	for kua_gene_boot,vua_gene_boot in ua_gene_boot_dict.iteritems():
		up_adj_boot.write('\t'.join([kua_gene_boot,str(vua_gene_boot)])+'\n')
	

	dn_adj_boot=open(str(args.geneid)+'_over_'+str(args.limit)+'_FPKM_in_'+str(args.no_cells)+'_Cells_'+str(args.iterations)+'iterations_dn_adj_bt.txt', 'w')
	for kda_gene_boot,vda_gene_boot in da_gene_boot_dict.iteritems():
		dn_adj_boot.write('\t'.join([kda_gene_boot,str(vda_gene_boot)])+'\n')
	
	#  p-values & logFC
	
	up_out_pvallogFC= open(str(args.geneid)+'_over_'+str(args.limit)+'_FPKM_in_'+str(args.no_cells)+'_Cells_'+str(args.iterations)+'iterations_up_pv.txt', 'w')
	for ku_gene_pvallogFC,vu_gene_pvallogFC in u_gene_pvallogFC_dict.iteritems():
		up_out_pvallogFC.write('\t'.join([str(seed), ku_gene_pvallogFC, '\t'.join(vu_gene_pvallogFC[0]) ])+'\n')
	
	
	dn_out_pvallogFC=open(str(args.geneid)+'_over_'+str(args.limit)+'_FPKM_in_'+str(args.no_cells)+'_Cells_'+str(args.iterations)+'iterations_dn_pv.txt', 'w')
	for kd_gene_pvallogFC, vd_gene_pvallogFC in d_gene_pvallogFC_dict.iteritems():
		dn_out_pvallogFC.write('\t'.join([str(seed), kd_gene_pvallogFC, '\t'.join(vd_gene_pvallogFC[0]) ])+'\n')
	
	up_adj_pvallogFC=open(str(args.geneid)+'_over_'+str(args.limit)+'_FPKM_in_'+str(args.no_cells)+'_Cells_'+str(args.iterations)+'iterations_up_adj_pv.txt', 'w')
	for kua_gene_pvallogFC, vua_gene_pvallogFC in ua_gene_pvallogFC_dict.iteritems():
		up_adj_pvallogFC.write('\t'.join([str(seed), kua_gene_pvallogFC, '\t'.join(vua_gene_pvallogFC[0]) ])+'\n')
	
	dn_adj_pvallogFC=open(str(args.geneid)+'_over_'+str(args.limit)+'_FPKM_in_'+str(args.no_cells)+'_Cells_'+str(args.iterations)+'iterations_dn_adj_pv.txt', 'w')
	for kda_gene_pvallogFC, vda_gene_pvallogFC in da_gene_pvallogFC_dict.iteritems():
		dn_adj_pvallogFC.write('\t'.join([str(seed), kda_gene_pvallogFC, '\t'.join(vda_gene_pvallogFC[0]) ])+'\n')
	