read.csv.list <- function(file,comment.char='#',...) {
    text <- readLines(file)
    mm <- vector('integer',length(text))
    mm[c(0,grep("^[[:space:]]*$",text))] <- 1
    mm <- cumsum(mm)
    
    tables <- lapply(split(text,mm),
        function (block) { read.csv(textConnection(block),comment.char=comment.char,...) })
    names(tables) <- lapply(split(text,mm),
        function (block) {
            fl <- scan(textConnection(block),n=1,blank.lines.skip=TRUE,sep='\n',what='');
            if (substr(fl,1,1)=='#') sub('^#[[:space:]]*','',fl) else ''; })
    
    return(tables)
}

