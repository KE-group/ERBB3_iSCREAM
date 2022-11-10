#annovar_cDNA_Find# Function Description

# This function takes a mutation info column from ANNOVAR output and selects the cDNA change based on the isoform provided.

# Isoform accepted as GI identifier in the block for the particular isoform in the ANNOVAR mutation info. You should look at the data and then provide the info. Defaults to ERBB4 A2 isoform.

annovar_cDNA_Find=function(MutationColumn,isoform="NM_001042599"){
  MutationList=c("List of mutations")
  for(i in seq(1:length(MutationColumn))){
    MutInfo=MutationColumn[i]
    l=sort(unique(unlist(strsplit(MutInfo,","))))
    l2=l[grep(isoform,l)]
    l2.s=unique(unlist(strsplit(l2,":")))
    l3=l2.s[grep("^c",l2.s)]
    l3=gsub("c.","",l3)
    #l3=gsub("X","*",l3)
    MUTATION=l3
    if(length(MUTATION)==0) MUTATION=" "
    MutationList=c(MutationList,MUTATION)
  }
  return(MutationList[-1])
}


# MutationColumn=Mutation.Table$AAChange.refgene;i=7591;isoform="NM_005228"
# annovarMutCodeFind(tab.s$Mutation)
