library(httr)
library(jsonlite)
library(RCurl)
httr::set_config(config(ssl_verifypeer = 0L))

#URL base for GeneCodis
urlBase = "https://genecodis.genyo.es/gc4/"
checkRatelimit <- function(content){
  ratelimit=FALSE
  #print("checking rate-limit...")
  if(length(grep("rate-limit",content))>0){
    ratelimit=TRUE
  }
  return(ratelimit)
}

#function to launch analysis
launchAnalysis <- function(organism,inputType,inputQuery,annotationsDBs,enrichmentStat="hypergeom",
                           universeScope="annotated",inputCoannotation="coannotation_no",coannotationAlgorithm="fpgrowth",minimalInputCoannotation=10,
                           secondInputQuery=list(),inputName1="input1",inputName2="input2",customUniverse=list(),inputEmail="",ReportName="input1"){
  
  params <- list(inputmode = "on",
                 email = inputEmail,
                 inputSupport = minimalInputCoannotation,
                 input = list(input=list(inputQuery),input2=list(secondInputQuery)),
                 jobName= ReportName,
                 annotations = list(annotationsDBs),
                 stat = enrichmentStat,
                 scope = universeScope,
                 algorithm = coannotationAlgorithm,
                 inputNames = list(input=inputName1, input2=inputName2),
                 inputtype = inputType,
                 organism = organism,
                 universe = customUniverse,
                 gc4uid = "",
                 coannotation = inputCoannotation)
  #convert list of params to json
  jsonparams <- jsonlite::toJSON(params,auto_unbox=T,pretty = F)
  #define URL to make request
  requestJobURL = paste0(urlBase,"analysis")
  #make request to GeneCodis
  print('Got analysis petition...')
  print('Performing the analyses...')
  request<-httr::POST(url = requestJobURL, body=jsonparams, encode="json", content_type("application/json"))
  #check if analysis has finished
  request_content=content(request,"text")
  if(status_code(request)!=200){
    stop(request_content)
  }
  params$gc4uid=substr(x = request_content,start = gregexpr(":", request_content)[[1]]+1,stop = nchar(request_content) )
  #check invalid input
  stateURL = paste0(urlBase,"analysisResults?job=",params$gc4uid)
  #make request to GeneCodis
  request = GET(stateURL)
  status=content(GET(paste0(urlBase,"checkstate?job=",params$gc4uid)))
  while(status$state=="PENDING"){
    status=content(GET(paste0(urlBase,"checkstate?job=",params$gc4uid)))
  }
  request = GET(stateURL)
  #check if error happens and stop execution
  warn_for_status(request)
  stop_for_status(request)
  #get content of request
  request_content=content(request,"text")
  #check if user reached rate limit
  if(checkRatelimit(request_content)==TRUE){
    stop("You have exceeded the GeneCodis rate-limit, which is 10 requests per minute. If you think you may need a further access to GeneCodis please contact us at bioinfo@genyo.es")
  }
  #check if input is invalid
  if(length(grep("invalid input list",request_content))>0){
    stop("It seems that you have provided an invalid input list. It could be also possible that the input does not match the selected organism or that we have no record of those in our database associated to the selected annotations. Sorry the inconveniences. Please go to the Help tab and check the allowed ids.")
  }
  #check if unexpected error happened
  if(length(grep("unexpected",request_content))>0){
    stop(paste0("Please send the job ticket, ",params$gc4uid," to bioinfo@genyo.es, the server found an unexpected error. We will solve it as soon as possible."))
  }
  print(paste0('Analysis finished, jobID:',params$gc4uid))
  print('Generating results...')
  #make request to obtain results
  return(getResults(params$gc4uid))
}

getResults <- function(gc4uid){
  stats_tables=list()
  results=list()
  results["jobID"]=gc4uid
  resultsurl=paste0(urlBase,"results?job=",gc4uid,"&annotation=all&jsonify=t")
  request=httr::GET(resultsurl)
  #check if error happens and stop execution
  warn_for_status(request)
  stop_for_status(request)
  result=content(request,"text")
  result <- fromJSON(result)
  for(lista in names(result)){
    temp1=t.data.frame(data.frame(matrix(unlist(result[[lista]]), nrow=length(result[[lista]]), byrow=TRUE)))
    colnames(temp1)=c("annotation_id","annotsbias","description","genes","genes_found","input_size","pval","pval_adj","relative_enrichment","term_genes","universe")
    rownames(temp1)=c(1:nrow(temp1))
    temp1=data.frame(temp1)
    stats_tables[[lista]]=temp1
  }
  results[["stats_tables"]]=stats_tables
  resultsurl=paste0(urlBase,"qc?job=",gc4uid)
  request=httr::GET(resultsurl)
  #check if error happens and stop execution
  warn_for_status(request)
  stop_for_status(request)
  result=content(request,"text")
  result=fromJSON(result)
  results[["quality_controls"]]=c(result)
  return(results)
}

queryGenes <- function(genes,orgs,databases){
  print("Quering genecodis database...")
  listgenes=""
  for(gene in genes){
    listgenes=paste0(listgenes,gene,",")
  }
  resultsurl=paste0(urlBase,"queryGene?genes=",listgenes)
  print(resultsurl)
  if(!missing(orgs)){
    listorgs=""
    for(org in orgs){
      listorgs=paste0(listorgs,org,",")
    }
    resultsurl=paste0(resultsurl,"&orgs=",listorgs)
  }
  print(resultsurl)
  if(!missing(databases)){
    listdbs=""
    for(db in databases){
      listdbs=paste0(listdbs,db,",")
    }
    resultsurl=paste0(resultsurl,"&databases=",listdbs)
    resultsurl=substring(resultsurl,1,nchar(resultsurl)-1)
  }
  print(resultsurl)
  request=httr::GET(resultsurl)
  #check if error happens and stop execution
  warn_for_status(request)
  stop_for_status(request)
  result=content(request,"parsed")
  table=read.csv(text=result,header = TRUE,sep="\t")  
  print("Query finished")
  return(table)
}


