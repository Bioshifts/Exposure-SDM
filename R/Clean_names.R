# Clean scientific names to increase probability of resolving names
require(stringr)
Clean_Names <- function(spname, return_gen_sps = FALSE){
    pbsapply(spname, .clean_Names, return_gen_sps = return_gen_sps)
}
.clean_Names <- function(spname_i, return_gen_sps = FALSE){
    if(is.na(spname_i) | spname_i == "NA"){
        return("NA")
    }
    else{
        # replace artifacts from encoding issue
        tmp. <- gsub('<U\\+FB01>', 'fi', spname_i)
        tmp. <- gsub('ÿ', '', tmp.)
        # dot within name
        tmp. <- strsplit(tmp.,"[.]")[[1]]
        tmp. = paste(tmp., collapse = " ")
        # reduces repeated white space
        tmp. <- str_squish(tmp.)
        # remove letters with accent
        tmp. <- iconv(tmp., "latin1", "ASCII", "")
        # remove parentheses
        tmp. <- strsplit(tmp., " ")[[1]]
        if(any(grepl("[(]", tmp.))){
            tmp. = tmp.[-which(grepl("[(]", tmp.))]
        }
        # remove bar
        if(any(grepl("[/]", tmp.))){
            pos <- grep("[/]", tmp.)
            tmp.[pos] = gsub("[/]"," ",tmp.[pos])
            tmp. = paste(tmp., collapse = " ")
            tmp. <- strsplit(tmp., " ")[[1]]
        }
        # Genus uppercase and species lowercase
        tmp.1 = str_to_title(tmp.[1])
        if(length(tmp.)>1){
            tmp.2 <- sapply(tmp.[-1], function(x) str_to_lower(x))
            tmp.2 = paste(tmp.1, paste(tmp.2, collapse = " "), collapse = " ")
            tmp. <- strsplit(tmp.2, " ")[[1]]
        } else {
            tmp. <- tmp.1
        }
        # fix sp. var. subspecies
        tmp.[which(tmp. == "sp")] <- "sp."
        tmp.[which(tmp. == "spp")] <- "sp."
        tmp.[which(tmp. == "spp.")] <- "sp."
        tmp.[which(tmp. == "var")] <- "var."
        tmp.[which(tmp. == "cf")] <- "cf."
        tmp.[which(tmp. == "ssp")] <- "subsp."
        tmp.[which(tmp. == "ssp.")] <- "subsp."
        tmp.[which(tmp. == "subsp")] <- "subsp."
        tmp. <- gsub("[..]",".",tmp.)
        if(return_gen_sps) {
            # remove var subsp
            if(any(tmp. == "var" | tmp. == "subsp" | tmp. == "var." | tmp. == "cf." | tmp. == "subsp." | tmp. == "X" | tmp. == "x" | tmp. == "[,]")){
                tmp. = tmp.[-which(tmp. == "var" | tmp. == "subsp" | tmp. == "var." | tmp. == "cf." | tmp. == "subsp." | tmp. == "X" | tmp. == "x" | tmp. == "[,]")]
            }
            if(length(tmp.)>2){
                tmp. <- tmp.[1:2]
            }
        }
        if(length(tmp.)==1){
            tmp. = paste(tmp., "sp.", collapse = " ")
        } else {
            tmp. = paste(tmp., collapse = " ")
        }
        
        return(tmp.)
    }
}