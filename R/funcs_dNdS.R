#' @title Function to calculate the synonymous vs nonsynonymous substitutionrate for two organisms.
#' @description This function takes 
#' @param query_file a character string specifying the path to the CDS file of interest (query organism).
#' @param subject_file a character string specifying the path to the CDS file of interest (subject organism).
#' @param blast_path a character string specifying the path to the BLAST program (in case you don't use the default path).
#' @param multialn_path a character string specifying the path to the multiple alignment program (in case you don't use the default path).
#' @param codonaln_path a character string specifying the path to the codon alignment program (in case you don't use the default path).
#' @param comp_cores a numeric value specifying the number of cores that shall be
#' @param tool a character string specifying the program that should be used e.g. "clustalw". 
#' @author Sarah Scharfenberg and Hajk-Georg Drost 
#' @details This function...
#' @return A data.table storing the dNdS values of the correspnding genes.
#' @import data.table
#' @export
dNdS <- function(query_file, subject_file, 
                 blast_mode="best hit", blast_path = NULL, 
                 multialn_tool="clustalw", multialn_path = NULL,
                 codonaln_tool="pal2nal", codonaln_path = NULL,
                 dnds_tool="gestimator", dnds_path = NULL, comp_cores = 1){
        
        if(!is.element(blast_mode, c("best hit","recursive")))
                stop("Please choose a blast mode that is supported by this function.")
        
        if(!is.element(multialn_tool, c("clustalw", "clustalo","muscle", "tcoffee")))
                stop("Please choose a multiple alignment tool that is supported by this function.")
        
        if(!is.element(codonaln_tool, c("pal2nal")))
                stop("Please choose a codon alignment tool that is supported by this function.")
        
        if(!is.element(dnds_tool, c("gestimator")))
                stop("Please choose a dnds tool that is supported by this function.")
        
        # blast each translated aminoacid sequence against the related database to get a 
        # hit table with pairs of geneids  
        
        query_cds <- NULL
        query_aa <- NULL
        subject_cds <- NULL
        subject_aa <- NULL
        
        if(blast_mode == "best hit"){
#                 
#                 hit.table <- data.table::copy(
#                         blast_best(query_file = query_file, subject_file = subject_file, 
#                                    path = blast_path, comp_cores = comp_cores))
#                 
                hit.table <- blast_best(query_file = query_file, subject_file = subject_file, 
                                                   path = blast_path, comp_cores = comp_cores)
                                
                query_cds <- read.cds(file = query_file, format = "fasta")
                subject_cds <- read.cds(file = subject_file, format = "fasta")
                
                query_aa <- read.proteome(file = "_blast/blastinput.fasta", format = "fasta")
                
                filename <- unlist(strsplit(subject_file, "/", fixed = FALSE, perl = TRUE, useBytes = FALSE))
                filename <- filename[length(filename)]
                subject_aa <- read.proteome(file = paste0("_database/out_",filename,"_translate.fasta"), format = "fasta")
    
        }
        
        if(blast_mode == "recursive"){
                               
                hit.table <- data.table::copy(
                        blast_rec(query_file = query_file, subject_file = subject_file, 
                                   path = blast_path, comp_cores = comp_cores))
                
                query_cds <- read.cds(file = query_file, format = "fasta")
                subject_cds <- read.cds(file = subject_file, format = "fasta")
                
                
                filename <- unlist(strsplit(query_file, "/", fixed = FALSE, perl = TRUE, useBytes = FALSE))
                filename <- filename[length(filename)]
                query_aa <- read.proteome(file = paste0("_database/out_",filename,"_translate.fasta"), format = "fasta")
                
                filename <- unlist(strsplit(subject_file, "/", fixed = FALSE, perl = TRUE, useBytes = FALSE))
                filename <- filename[length(filename)]
                subject_aa <- read.proteome(file = paste0("_database/out_",filename,"_translate.fasta"), format = "fasta")
                
                
        }
     
        data.table::setnames(query_cds, old=c("geneids", "seqs"), new = c("query_id","query_cds"))
        data.table::setnames(subject_cds, old=c("geneids", "seqs"), new = c("subject_id","subject_cds"))
        data.table::setnames(query_aa, old=c("geneids", "seqs"), new = c("query_id","query_aa"))
        data.table::setnames(subject_aa, old=c("geneids", "seqs"), new = c("subject_id","subject_aa"))             
        
        hit.table <- dplyr::inner_join(hit.table, query_cds, by = "query_id")
        hit.table <- dplyr::inner_join(hit.table, query_aa, by = "query_id")
        hit.table <- dplyr::inner_join(hit.table, subject_cds, by = "subject_id")
        hit.table <- dplyr::inner_join(hit.table, subject_aa, by = "subject_id")
      
        
        # foreach hit-pair do
        # - create a fasta file containing the aminoacid sequences -> aa.fasta
        # - create a fatsa file containing the cds -> cds.fasta
        # - do multi_aln on aa.fasta -> aa.aln
        # - do codon_aln on aa.aln and cds.fasta -> codon.aln
        # - do compute_dnds on codon.aln -> fill in hittable
        
      #hit.table[, dnds:= compute_dnds( <Rest der Spalte ?? >), by=c("query_id", "subject_id")]
        
      # obacht hit.table enthalt noch den evalue
      hit.table <- data.table(hit.table)
      data.table::setkey(hit.table, query_id)
      hit.table[, dnds:=as.vector(apply(.SD,1, FUN=function(x){ compute_dnds(x,
                               multialn_tool = "clustalw", codonaln_tool = "pal2nal", 
                               dnds_tool = "gestimator",
                               codonaln_path = "/home/sarah/Programs/pal2nal.v14/" )})),]
      return(hit.table)  
  
      
   #   compute_dnds(res[1,], 
   #                multialn_tool = "clustalw", codonaln_tool = "pal2nal", 
   #                dnds_tool = "gestimator",
   #                codonaln_path = "/home/sarah/Programs/pal2nal.v14/")

    #  res[,apply(.SD,1,function(x) View(x) ), by=query_id]
    # by query id results in a frame with id and lists of 6
}

#' @title Function to calculate the synonymous vs nonsynonymous substitutionrate for a codon alignment.
#' @description This function takes 
#' @param file a character string specifying the path to a codon alignment file 
#' @author Sarah Scharfenberg and Hajk-Georg Drost 
#' @details This function ...
#' @return some value
#' @import data.table
#' @export
substitutionrate <- function(file, tool, path = NULL){
        
        if(!is.element(tool,c("gestimator")))
                stop("Please choose a tool that is supported by this function.")
        
        if(!file.exists("_calculation/")){
                
                dir.create("_calculation")
        }
        
        if(tool == "gestimator"){
                
                # file in fasta required
                
            tryCatch(
            {    
                system(paste0("gestimator -i ",file," -o _calculation/gestimout"))
                hit.table <-data.table::fread("_calculation/gestimout")
                data.table::setnames(hit.table, old=c("V1","V2","V3","V4","V5"), 
                                     new = c("query_id","subject_id","dN","dS","dNdS"))
                data.table::setkey(hit.table, query_id)
                return(hit.table)
            },error = function(){ print(paste0("Please check the correct path to ",tool,
                                               "... the interface call did not work properly.") ) }
            
            )
        }
        
}

compute_dnds <- function(x, 
                         multialn_tool="clustalw", multialn_path = NULL,
                         codonaln_tool="pal2nal", codonaln_path = NULL,
                         dnds_tool="gestimator", dnds_path = NULL){
 
        names <- list(x[[1]],x[[2]])
        seqs <- list(x[[4]],x[[6]])
        aa <- list(x[[5]],x[[7]])

        # create cds fasta
        seqinr::write.fasta(sequences = seqs, names = names, file.out = "_alignment/cds.fasta")
        
        # create aa fasta
        seqinr::write.fasta(sequences = aa, names = names, file.out = "_alignment/aa.fasta")
        
        # align aa -> <multialn_tool>.aln
        multi_aln(file = "_alignment/aa.fasta", 
                  tool = multialn_tool, 
                  get_aln = FALSE, path = multialn_path)

        # align codon -> cds.aln
        codon_aln(file_aln = paste0("_alignment/",multialn_tool,".aln"),
                  file_nuc = "_alignment/cds.fasta",
                  tool = codonaln_tool,
                  format = "fasta",
                  get_aln = FALSE, path = codonaln_path)
        
        # compute kaks
        hit.table <- substitutionrate(file = paste0("_alignment/",codonaln_tool,".aln"), 
                                      tool = "gestimator", path = dnds_path)
        return(hit.table[,dNdS])

}