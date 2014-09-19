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
        
        if(!is.element(dnds_tool, c("gestimator","li")))
                stop("Please choose a dnds tool that is supported by this function.")
        
        # blast each translated aminoacid sequence against the related database to get a 
        # hit table with pairs of geneids  
        
  
        if(blast_mode == "best hit"){
                
                hit.table <- data.table::copy(
                        blast_best(query_file = query_file, subject_file = subject_file, 
                                   path = blast_path, comp_cores = comp_cores))
                
#                 hit.table <- blast_best(query_file = query_file, subject_file = subject_file, 
#                                                    path = blast_path, comp_cores = comp_cores)
#                                 
                q_cds <- read.cds(file = query_file, format = "fasta")
                s_cds <- read.cds(file = subject_file, format = "fasta")
                
                q_aa <- read.proteome(file = "_blast/blastinput.fasta", format = "fasta")
                
                filename <- unlist(strsplit(subject_file, "/", fixed = FALSE, perl = TRUE, useBytes = FALSE))
                filename <- filename[length(filename)]
                s_aa <- read.proteome(file = paste0("_database/out_",filename,"_translate.fasta"), format = "fasta")
    
        }
        
        if(blast_mode == "recursive"){
                               
                hit.table <- data.table::copy(
                        blast_rec(query_file = query_file, subject_file = subject_file, 
                                   path = blast_path, comp_cores = comp_cores))
                
                q_cds <- read.cds(file = query_file, format = "fasta")
                s_cds <- read.cds(file = subject_file, format = "fasta")
                
                
                filename <- unlist(strsplit(query_file, "/", fixed = FALSE, perl = TRUE, useBytes = FALSE))
                filename <- filename[length(filename)]
                q_aa <- read.proteome(file = paste0("_database/out_",filename,"_translate.fasta"), format = "fasta")
                
                filename <- unlist(strsplit(subject_file, "/", fixed = FALSE, perl = TRUE, useBytes = FALSE))
                filename <- filename[length(filename)]
                s_aa <- read.proteome(file = paste0("_database/out_",filename,"_translate.fasta"), format = "fasta")
                
                
        }
    

        data.table::setnames(q_cds, old=c("geneids", "seqs"), new = c("query_id","query_cds"))
        data.table::setnames(s_cds, old=c("geneids", "seqs"), new = c("subject_id","subject_cds"))
        data.table::setnames(q_aa, old=c("geneids", "seqs"), new = c("query_id","query_aa"))
        data.table::setnames(s_aa, old=c("geneids", "seqs"), new = c("subject_id","subject_aa"))
        
#### A #####

        hit.table <- dplyr::inner_join(tbl_dt(hit.table), tbl_dt(s_cds), by = "subject_id")
        hit.table <- dplyr::inner_join(tbl_dt(hit.table), tbl_dt(s_aa), by = "subject_id")
        hit.table <- dplyr::inner_join(tbl_dt(hit.table), tbl_dt(q_cds), by = "query_id")
        hit.table <- dplyr::inner_join(tbl_dt(hit.table), tbl_dt(q_aa), by = "query_id")  
        
#           tbl_cds <- dplyr::inner_join(dplyr::tbl_dt(q_cds), dplyr::tbl_dt(s_cds), by = c("query_id","subject_id"))
#           tbl_aa <- dplyr::inner_join(dplyr::tbl_dt(q_aa), dplyr::tbl_dt(s_aa), by = c("query_id","subject_id"))
#           tbl_cds_aa <- dplyr::inner_join(dplyr::tbl_dt(tbl_cds), dplyr::tbl_dt(tbl_aa), by = "query_id"))
#           final_tbl <- dplyr::inner_join(dplyr::tbl_dt(hit.table), dplyr::tbl_dt(tbl_cds_aa), by = "query_id"))
#           hit.table <- data.table::copy(final_tbl)

# AN DIESER STELLE MUSS ICH MICH DOCH SEHR WUNDERN     
# DURCH DIE ÜBERGABE DES INNER ERGEBNISSES DES INNER JOIN IST DAS RESULTAT
# HIER EIN DATA.TABLE DER PLÖTZLICH IN DIE KEY SPALTEN GETAUSCHT HAT!
# ALLERDINGS IST DAS NICHT MIT SEKEY() RÜCKGÄNGIG ZU MACHEN
# AUẞERDEM IST DER KEY NUR NOCH SUBJECT_ID



### B ###

# Variante B mit hit.table[q_cds] läuft am Ende auch auf das Problem hinaus, dass setkey
# für 2 Aufrufe umbesetzt werden muss.
# Auch muss dabei ebenfalls hit.table überschrieben werden    

        #data.table::setkeyv(hit.table,c("query_id", "subject_id"))

# foreach hit-pair do
        # - create a fasta file containing the aminoacid sequences -> aa.fasta
        # - create a fatsa file containing the cds -> cds.fasta
        # - do multi_aln on aa.fasta -> aa.aln
        # - do codon_aln on aa.aln and cds.fasta -> codon.aln
        # - do compute_dnds on codon.aln -> fill in hittable
        
      #hit.table[, dnds:= compute_dnds( <Rest der Spalte ?? >), by=c("query_id", "subject_id")]
        
      # obacht hit.table enthalt noch den evalue


      hit.table[, dnds:=as.vector(apply(.SD, 1 ,FUN=function(x){ compute_dnds(x,
                               multialn_tool = multialn_tool, codonaln_tool = codonaln_tool, 
                               dnds_tool = dnds_tool,
                               codonaln_path = codonaln_path )}))]

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
#' @param format a character string specifying the format a used in read.alignment (for Li) 
#' @author Sarah Scharfenberg and Hajk-Georg Drost 
#' @details This function ...
#' @return some value
#' @import data.table
#' @export
substitutionrate <- function(file, tool, format="fasta", path = NULL){
        
        if(!is.element(tool,c("gestimator","Li")))
                stop("Please choose a tool that is supported by this function.")
        
        if(!is.element(format,c("mase", "clustal", "phylip", "fasta" , "msf" )))
                stop("Please choose a format that is supported by seqinr::read.alignment.")
        
        if(!file.exists("_calculation/")){
                
                dir.create("_calculation")
        }
        
        if(tool == "gestimator"){
                
                # file in fasta required
                if(format != "fasta")
                        stop("To use gestimator fasta format is required.")
                
            tryCatch(
            {    
                system(paste0("gestimator -i ",file," -o _calculation/gestimout"))
                hit.table <-data.table::fread("_calculation/gestimout")
                data.table::setnames(hit.table, old=c("V1","V2","V3","V4","V5"), 
                                     new = c("query_id","subject_id","dN","dS","dNdS"))
                data.table::setkey(hit.table, query_id)
                
                print("Substitutionrate successfully calculated by gestimator")
                
                return(hit.table)
            },error = function(){ print(paste0("Please check the correct path to ",tool,
                                               "... the interface call did not work properly.") ) }
            
            )
        }
        
        if(tool == "Li" ){
                
                        aln <- seqinr::read.alignment(file = file, format = format)
                
                        res_list <- seqinr::kaks(aln)
                        
                        res <- data.table(t(c(unlist(aln$nam),  res_list$ka[1], res_list$ks[1], (res_list$ka[1]/res_list$ks[1]))))
                        data.table::setnames(res, old = paste0("V",1:5), 
                                             new = c("query_id","subject_id", "dN", "dS","dNdS"))
                 
                        print("Substitutionrate successfully calculated using Li's method.")
                        
                        return(res)
                        
        
        }
        
}

compute_dnds <- function(x, 
                         multialn_tool="clustalw", multialn_path = NULL,
                         codonaln_tool="pal2nal", codonaln_path = NULL,
                         dnds_tool="gestimator", dnds_path = NULL){
        #return(x["query_id"])
        names <- list(x["query_id"],x["subject_id"])
        seqs <- list(x["query_cds"],x["subject_cds"])
        aa <- list(x["query_aa"],x["subject_aa"])

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
                  get_aln = FALSE)
        
        # compute kaks
        hit.table <- substitutionrate(file = paste0("_alignment/",codonaln_tool,".aln"), 
                                      tool = dnds_tool, path = dnds_path)
        return(hit.table[,dNdS])

}
# 
# When running
# table <- dNdS("data/ortho_thal_cds.fasta", "data/ortho_lyra_cds.fasta")
# I still get this warning and dont know how to fix.
# 
# Warning message:
# In `[.data.table`(hit.table, , `:=`(dnds, as.vector(apply(.SD, 1,  :
# Invalid .internal.selfref detected and fixed by taking a copy of the whole 
# table so that := can add this new column by reference. At an earlier point, 
# this data.table has been copied by R (or been created manually using structure() 
# or similar). Avoid key<-, names<- and attr<- which in R currently (and oddly) 
# may copy the whole data.table. Use set* syntax instead to avoid copying: ?set, 
# ?setnames and ?setattr. Also, in R<=v3.0.2, list(DT1,DT2) copied the entire DT1 
# and DT2 (R's list() used to copy named objects); please upgrade to R>v3.0.2 
#  if that is biting. If this message doesn't help, please report to 
#  datatable-help so the root cause can be fixed.
