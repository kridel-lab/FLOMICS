#####################################################################################################
#functions
#1 Dividing datasets based on their Chromosome number
array2list <-function(data, column='subtype') {
  newlist = lapply(unique(data[,column]), function(x) {
    A = data[which(data[,column]==x),, drop=FALSE]
    #rownames(A) = rownames(data)[data[,column]==x]
    return(A)
  })
  names(newlist) = unique(data[,column])
  return(newlist)
}

#####################################################################################################
#1.2 sort a sequence of characters
sortection=function(seq.list){
  seq.list.new=lapply(seq.list, function(x){
    df=data.frame(seq=x,ord=1:length(x))
    return(df)
  })
  #Take the intersections
  Union <- Reduce(function(x, y) merge(x, y, all=FALSE, by = "seq"), seq.list.new, accumulate=F)
  
  #prepare the output df
  df=Union[,-1]
  rownames(df)=Union[,1]
  
  return(df)
}

#####################################################################################################
Bin.No=function(x,a1,a2,a3,a4){
  if(x<=a1) {
    return(1)
  } else if(x>a1&x<=a2) {
    return(2)
  } else if(x>a2&x<=a3) {
    return(3)
  } else if(x>a3&x<=a4){
    return(4)
  }
}
