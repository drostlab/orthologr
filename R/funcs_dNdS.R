dNdS <- function(query_file, subject_file, blast_mode="best hit", multialn_tool="clustalw", 
                 codonaln_tool="pal2nal", dnds_tool="gestimator", path = NULL, comp_cores = 1){
        
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
                
                hit.table <- data.table::copy(
                        blast_best(query_file = query_file, subject_file = subject_file, 
                                   path = path, comp_cores = comp_cores))
                
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
                                   path = path, comp_cores = comp_cores))
                
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
        
        hit.table <- dplyr::inner_join( hit.table, query_cds, by = "query_id", )
        hit.table <- dplyr::inner_join(hit.table, query_aa, by = "query_id")
        hit.table <- dplyr::inner_join(hit.table, subject_cds, by = "subject_id")
        hit.table <- dplyr::inner_join(hit.table, subject_aa, by = "subject_id")
      
      return(hit.table)  
        
        # foreach hit-pair do
        # - create a fasta file containing the aminoacid sequences -> aa.fasta
        # - create a fatsa file containing the cds -> cds.fasta
        # - do multi_aln on aa.fasta -> aa.aln
        # - do codon_aln on aa.aln and cds.fasta -> codon.aln
        # - do compute_dnds on codon.aln -> fill in hittable
        
      #hit.table[, dnds:= compute_dnds(codon_aln(multi_aln(aa.fasta), cds.fasta)), by=c("query_id", "subject_id")]
        
        

}