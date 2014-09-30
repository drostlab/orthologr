

clean_all_folders <- function(){
        
        if(file.exists("_alignment"))
                file.remove("_alignment")
        
        if(file.exists("_blast_db"))
                file.remove("_blast_db")
        
        if(file.exists("_calculation"))
                file.remove("_calculation")
                
}