#### Phylogenetic PATH analysis
#### Follows approach of Ch 8 in Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology

require(ape)
require(caper)
require(nlme)
require(geiger)


###setting up the data Matrix###

war_PA<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Evolution\ figures/Parulidae_PA.csv")

war_PA[is.na(war_PA)]<-0

sex_Dim<-read.csv("~/Dropbox/Warbler.Molt.Migration/Parulidae_sexual_dimorphism.csv")
sex_Dim[is.na(sex_Dim)]<-0
sextent<-sex_Dim[,2:13];rownames(sextent)<-sex_Dim[,1]
sDim<-rowSums(sextent);names(sDim)<-rownames(sextent)


temp_Dim<-read.csv("~/Dropbox/Warbler.Molt.Migration/Parulidae_temporal_dimorphism.csv")
temp_Dim[is.na(temp_Dim)]<-0
textent<-temp_Dim[,2:13];rownames(textent)<-temp_Dim[,1]
tDim<-rowSums(textent);names(tDim)<-rownames(textent)

war.mig<-read.table("~/Dropbox/Warbler.Molt.Migration/migratory_distance.txt")


PA<-data.frame(war_PA$tree_taxon,war_PA$extent);names(PA)<-c("Birdlife","PA")
sexDim<-data.frame(names(sDim),rowSums(sextent));names(sexDim)<-c("Birdlife","sexDim")
seaDim<-data.frame(names(tDim),rowSums(textent));names(seaDim)<-c("Birdlife","seaDim")
mig<-data.frame(rownames(war.mig),war.mig$migdist1);names(mig)<-c("Birdlife","migDist")

p1<-merge(PA,sexDim,by="Birdlife")
p2<-merge(p1,seaDim,by="Birdlife")
p3<-merge(p2,mig,by="Birdlife")


source('~/Dropbox/Warbler.Molt.Migration/3_processing.R', chdir = TRUE)


###migdist comes from Processing.R
sol<-data.frame(migdist$AOU,migdist$daylength);names(sol)<-c("AOU","sol")




###--------------------------------------------------------------------------------#######
####Set all the taxon names to AOU names###
taxa.key=read.delim('~/Dropbox/Warbler.Molt.Migration/Parulidae_taxonomies.txt',stringsAsFactors=F)


name.check(p3$taxon,taxa.key$Birdlife)

p4<-merge(p3,taxa.key,by="Birdlife")
p5<-merge(p4,sol,by="AOU")

warData<-data.frame(p5$AOU,p5$PA,p5$sexDim,p5$seaDim,p5$migDist,p5$sol);names(warData)<-c("AOU","PA","SD","TD","Mig","sol")

war8<-read.nexus("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Parulidae_phylogeny/Trees/r8stree.tre")



####-----------------------------------------------------------------------------##########
###drop taxa not in the data set from tree
t<-PA$PA;names(t)<-p5$AOU
warTree<-treedata(war8,t)$phy


###Create caper Comparative data object####
com.dat<-comparative.data(warTree,p5,AOU,vcv=TRUE,vcv.dim=3,warn.dropped=TRUE)

###Outline conditional dependencies to make sure we are clear about hypothesis - five edges and eight vertices should produce six CDs

##function to get number of conditional dependencies by number of vertices (V) and edges (A)
condNum <- function(V, A) {
    (factorial(V)/(2 * factorial(V - 2))) - A
}


##Notes on translating directed acyclic graphs to conditional independencies: 
#first step is to use condNum to check how many CIs there should be
#second step is to list each unique set of non-adjacent vertices
#third step is to list parent variable of EACH non-adjecent variable - if variable has no parent don't list

##Note on translating CIs into Linear models
#start with the downstream non-adjacent vertex
# on the other side of the ~, list all other vertices in the CI, 
#last vertex should be the first vertex in the CI!! (this needs to be changed - also make sure "duplicated" equations are truly equivalent 


#Conditional independencies *indicates duplicated independencies
#Thes CI statements are based on the PATH hypotheses figure outlining the acyclic graph hypotheses

#Model 1 (6 CIs)

#(mig,pa){dl} 		pa~dl+mig
#(mig,td){pa}		td~pa+mig
#(mig,sd){td}		sd~td+mig
#(dl,td){pa,mig}	td~pa+mig+dl
#(dl,sd){td,mig}	sd~td+mig+dl
#(pa,sd){td,dl}		sd~td+dl+pa



#Model 2 (4 CIs) 

#*(mig,pa){dl} 		*pa~dl+mig
#*(mig,td){pa}		*td~pa+mig
#*(dl,sd){td,mig}	*sd~td+mig+dl
#*(pa,sd){td,dl}	*sd~td+dl+pa



#Model 3 (5 CIs)

#(mig,pa){td} 			pa~td+mig
#(mig,td){sd}			td~sd+mig
#(mig,dl){null set}		dl~mig
#(dl,pa){td}			pa~td+dl
#(sd,pa){td,mig}		pa~td+mig+sd


#Model 4 (5 CIs)

#*(mig,td){pa}		*td~pa+mig
#*(mig,sd){td}		*sd~td+mig
#*(dl,td){pa,mig}	*td~pa+mig+dl
#*(dl,sd){td,mig}	*sd~td+mig+dl
#*(pa,sd){td,dl}	*sd~td+dl+pa


#Model 5 (5 CIs)

#*(mig,pa){dl} 		*pa~dl+mig
#*(mig,td){pa}		*td~pa+mig
#*(mig,sd){td}		*sd~td+mig
#*(dl,sd){td,mig}	*sd~td+mig+dl
#*(pa,sd){td,dl}	*sd~td+dl+pa

#Model 6 (6 CIs)

#(mig,pa){td} 		pa~td+mig
#(mig,td){sd}		td~sd+mig
#(mig,sd){dl}		sd~dl+mig
#(dl,td){sd,mig}	td~sd+mig+dl
#(dl,pa){td,mig}	pa~td+mig+dl
#(sd,pa){td,dl}		pa~td+dl+sd


#Model 7 (6 CIs)

#*(mig,pa){dl} 		*pa~dl+mig
#(mig,td){dl}		td~dl+mig
#*(mig,sd){td}		*sd~td+mig
#(dl,sd){td}		sd~td+dl
#(pa,sd){td,dl}		sd~td+dl+pa
#(pa,td){dl}		td~dl+pa

#Model 8 (6 CIs)

#*(mig,pa){dl} 		*pa~dl+mig
#*(mig,td){pa}		*td~pa+mig
#*(dl,td){pa,mig}	td~pa+mig+dl
#(dl,sd){mig}		sd~mig+dl
#(pa,sd){mig,dl}	sd~mig+dl+pa
#(td,sd){pa,mig}	sd~pa+mig+td


#Model 9 (5 CIs)

#(mig,pa){td} 		pa~td+mig
#(mig,td){sd}		sd~td+mig
#(dl,td){sd,mig}	td~sd+mig+dl
#(dl,pa){mig}		sd~mig+dl
#(sd,pa){mig,td}	pa~td+mig+sd


#Model 10 (5 CIs)
#*(mig,td){pa,dl}	*td~pa+dl+mig
#*(mig,sd){td,dl}	*sd~td+dl+mig
#*(dl,sd){td}		*sd~td+dl
#*(dl,pa){mig}		*pa~mig+dl
#(pa,sd){td}		sd~td+pa



#Model 11 (5 CIs)
#*(dl,sd){mig}			*sd~mig+dl
#(dl,td){sd,pa}			td~sd+pa+dl
#(dl,pa){sd}			pa~sd+dl
#(mig,td){dl,sd,pa}		td~dl+sd+pa+mig
#(mig,pa){dl,sd}		pa~dl+sd+mig


 
#Model 12 (3 CIs)
#(pa){mig}		pa~mig
#(td){mig}		td~mig
#(sd){mig}		sd~mig

#Model 13 (4 CIs)

#*(mig,td){dl,pa}		*td~dl+pa+mig 
#*(mig,sd){dl,td} 		*sd~dl+td+mig
#*(dl,sd){td}			*sd~td+dl
#*(pa,sd){dl,td}		*sd~dl+td+pa

#Model 14 (5 CIs)
#*(mig,pa){dl}			*pa~dl+mig
#(mig,td){dl,pa}		td~dl+pa+mig
#(mig,sd){dl}			sd~dl+mig
#(dl,sd)				sd~dl
#(pa,sd){dl}			sd~dl+pa



#Model 15 (6 CIs)

#*(mig,pa){dl} 		*pa~dl+mig
#*(mig,td){pa}		*td~pa+mig
#*(dl,td){pa}		td~pa+dl
#(dl,sd){mig}		sd~mig+dl
#(pa,sd){mig,dl}	sd~mig+dl+pa
#(td,sd){pa,mig}	sd~pa+mig+td



###translate these conditional independencies into PGLMs

CICc <- function(C, q, n) {
    C + 2 * q * (n/(n - 1 - q))
}

##Coefficients vertex is the upstream non-adjacent vertex

#Model 1




m1.1<-pgls(PA~sol+migDist,data=com.dat,lambda="ML");m1.1p<-summary(m1.1)$coefficients["migDist",4]
m1.2<-pgls(seaDim~PA+migDist,data=com.dat,lambda="ML");m1.2p<-summary(m1.2)$coefficients["migDist",4]
m1.3<-pgls(sexDim~seaDim+migDist,data=com.dat,lambda="ML");m1.3p<-summary(m1.3)$coefficients["migDist",4]
m1.4<-pgls(seaDim~PA+migDist+sol,data=com.dat,lambda="ML");m1.4p<-summary(m1.4)$coefficients["sol",4]
m1.5<-pgls(sexDim~seaDim+migDist+sol,data=com.dat,lambda="ML");m1.5p<-summary(m1.5)$coefficients["sol",4]
m1.6<-pgls(sexDim~seaDim+sol+PA,data=com.dat,lambda="ML");m1.6p<-summary(m1.6)$coefficients["PA",4]

C1<--2*(log(m1.1p)+log(m1.2p)+log(m1.3p)+log(m1.4p)+log(m1.5p)+log(m1.6p))
C1.pval<-1-pchisq(C1,2*6)
CICc1<-CICc(C1,9,108)



#Model 2
m2.1<-m1.1
m2.2<-m1.2
m2.3<-m1.5
m3.4<-m1.6




C2<--2*(log(m1.1p)+log(m1.2p)+log(m1.5p)+log(m1.6p))
C2.pval<-1-pchisq(C2,2*4)
CICc2<-CICc(C2,11,108)



#Model 3



m3.1<-pgls(PA~seaDim+migDist,data=com.dat,lambda="ML");m3.1p<-summary(m3.1)$coefficients["migDist",4]
m3.2<-pgls(seaDim~sexDim+migDist,data=com.dat,lambda="ML");m3.2p<-summary(m3.2)$coefficients["migDist",4]
m3.3<-pgls(sol~migDist,data=com.dat,lambda="ML");m3.3p<-summary(m3.3)$coefficients["migDist",4]
m3.4<-pgls(PA~seaDim+sol,data=com.dat,lambda="ML");m3.4p<-summary(m3.4)$coefficients["sol",4]
m3.5<-pgls(PA~seaDim+migDist+sexDim,data=com.dat,lambda="ML");m3.5p<-summary(m3.5)$coefficients["sexDim",4]

C3<--2*(log(m3.1p)+log(m3.2p)+log(0.0001)+log(m3.4p)+log(m3.5p))
C3.pval<-1-pchisq(C3,2*5)
CICc3<-CICc(C3,10,108)




#Model 4



m4.1<-pgls(seaDim~PA+migDist,data=com.dat,lambda="ML");m4.1p<-summary(m4.1)$coefficients["migDist",4]
m4.2<-pgls(sexDim~seaDim+migDist,data=com.dat,lambda="ML");m4.2p<-summary(m4.2)$coefficients["migDist",4]
m4.3<-pgls(seaDim~PA+migDist+sol,data=com.dat,lambda="ML");m4.3p<-summary(m4.3)$coefficients["sol",4]
m4.4<-pgls(sexDim~seaDim+migDist+sol,data=com.dat,lambda="ML");m4.4p<-summary(m4.4)$coefficients["sol",4]
m4.5<-pgls(sexDim~seaDim+seaDim+sol+PA,data=com.dat,lambda="ML");m4.5p<-summary(m4.5)$coefficients["PA",4]

C4<--2*(log(m4.1p)+log(m4.2p)+log(m4.3p)+log(m4.4p)+log(m4.5p))
C4.pval<-1-pchisq(C4,2*5)
CICc4<-CICc(C4,10,108)


#Model 5


m5.1<-pgls(PA~sol+migDist,data=com.dat,lambda="ML");m5.1p<-summary(m5.1)$coefficients["migDist",4]
m5.2<-pgls(seaDim~PA+migDist,data=com.dat,lambda="ML");m5.2p<-summary(m5.2)$coefficients["migDist",4]
m5.3<-pgls(seaDim~PA+migDist,data=com.dat,lambda="ML");m5.3p<-summary(m5.3)$coefficients["migDist",4]
m5.4<-pgls(sexDim~seaDim+migDist+sol,data=com.dat,lambda="ML");m5.4p<-summary(m5.4)$coefficients["sol",4]
m5.5<-pgls(sexDim~seaDim+sol+PA,data=com.dat,lambda="ML");m5.5p<-summary(m5.5)$coefficients["PA",4]


C5<--2*(log(m5.1p)+log(m5.2p)+log(m5.3p)+log(m5.4p)+log(m5.5p))
C5.pval<-1-pchisq(C5,2*5)
CICc5<-CICc(C5,11,108)

#Model 6



m6.1<-pgls(PA~seaDim+migDist,data=com.dat,lambda="ML");m6.1p<-summary(m6.1)$coefficients["migDist",4]
m6.2<-pgls(seaDim~sexDim+migDist,data=com.dat,lambda="ML");m6.2p<-summary(m6.2)$coefficients["migDist",4]
m6.3<-pgls(sexDim~sol+migDist,data=com.dat,lambda="ML");m6.3p<-summary(m6.3)$coefficients["migDist",4]
m6.4<-pgls(seaDim~sexDim+migDist+sol,data=com.dat,lambda="ML");m6.4p<-summary(m6.4)$coefficients["sol",4]
m6.5<-pgls(PA~seaDim+migDist+sol,data=com.dat,lambda="ML");m6.5p<-summary(m6.5)$coefficients["sol",4]
m6.6<-pgls(PA~seaDim+sol+sexDim,data=com.dat,lambda="ML");m6.6p<-summary(m6.6)$coefficients["sexDim",4]

C6<--2*(log(m6.1p)+log(m6.2p)+log(m6.3p)+log(m6.4p)+log(m6.5p)+log(m6.6p))
C6.pval<-1-pchisq(C6,2*6)
CICc6<-CICc(C6,9,108)


#Model 7


m7.1<-m1.1
m7.2<-pgls(seaDim~sol+migDist,data=com.dat,lambda="ML");m7.2p<-summary(m7.2)$coefficients["migDist",4]
m7.3<-m1.3	
m7.4<-pgls(sexDim~seaDim+sol,data=com.dat,lambda="ML");m7.4p<-summary(m7.4)$coefficients["sol",4]
m7.5<-pgls(sexDim~seaDim+sol+PA,data=com.dat,lambda="ML");m7.5p<-summary(m7.5)$coefficients["PA",4]		
m7.6<-pgls(seaDim~sol+PA,data=com.dat,lambda="ML");m7.6p<-summary(m7.6)$coefficients["PA",4]	
			
C7<--2*(log(m1.1p)+log(m7.2p)+log(m1.3p)+log(m7.4p)+log(m7.5p)+log(m7.6p))
C7.pval<-1-pchisq(C7,2*6)
CICc7<-CICc(C7,9,108)

#Model 8



m8.1<-m1.1
m8.2<-m1.2
m8.3<-m1.4
m8.4<-pgls(sexDim~migDist+sol,data=com.dat,lambda="ML");m8.4p<-summary(m8.4)$coefficients["sol",4]
m8.5<-pgls(sexDim~migDist+sol+PA,data=com.dat,lambda="ML");m8.5p<-summary(m8.5)$coefficients["PA",4]		
m8.6<-pgls(sexDim~PA+migDist+seaDim,data=com.dat,lambda="ML");m8.6p<-summary(m8.6)$coefficients["seaDim",4]

C8<--2*(log(m1.1p)+log(m1.2p)+log(m1.4p)+log(m8.4p)+log(m8.5p)+log(m8.6p))
C8.pval<-1-pchisq(C8,2*6)
CICc8<-CICc(C8,9,108)


#Model 9




m9.1<-pgls(PA~seaDim+migDist,data=com.dat,lambda="ML");m9.1p<-summary(m9.1)$coefficients["migDist",4]
m9.2<-pgls(sexDim~seaDim+migDist,data=com.dat,lambda="ML");m9.2p<-summary(m9.2)$coefficients["migDist",4]
m9.3<-pgls(seaDim~sexDim+migDist+sol,data=com.dat,lambda="ML");m9.3p<-summary(m9.3)$coefficients["sol",4]
m9.4<-pgls(sexDim~migDist+sol,data=com.dat,lambda="ML");m9.4p<-summary(m9.4)$coefficients["sol",4]
m9.5<-pgls(PA~seaDim+migDist+sexDim,data=com.dat,lambda="ML");m9.5p<-summary(m9.5)$coefficients["sexDim",4]



C9<--2*(log(m9.1p)+log(m9.2p)+log(m9.3p)+log(m9.4p)+log(m9.5p))
C9.pval<-1-pchisq(C9,2*5)
CICc9<-CICc(C9,11,108)

#Model 10 



m10.1<-pgls(PA~seaDim+migDist,data=com.dat,lambda="ML");m10.1p<-summary(m10.1)$coefficients["migDist",4]
m10.2<-pgls(sexDim~seaDim+migDist,data=com.dat,lambda="ML");m10.2p<-summary(m10.2)$coefficients["migDist",4]
m10.3<-pgls(seaDim~sexDim+migDist+sol,data=com.dat,lambda="ML");m10.3p<-summary(m10.3)$coefficients["sol",4]
m10.4<-pgls(sexDim~migDist+sol,data=com.dat,lambda="ML");m10.4p<-summary(m10.4)$coefficients["sol",4]
m10.5<-pgls(PA~seaDim+migDist+sexDim,data=com.dat,lambda="ML");m10.5p<-summary(m10.5)$coefficients["sexDim",4]

C10<--2*(log(m10.1p)+log(m10.2p)+log(m10.3p)+log(m10.4p)+log(m10.5p))
C10.pval<-1-pchisq(C10,2*5)
CICc10<-CICc(C10,9,108)




#Model 11 (5 CIs)



m11.1<-pgls(sexDim~migDist+sol,data=com.dat,lambda="ML");m11.1p<-summary(m11.1)$coefficients["sol",4]
m11.2<-pgls(seaDim~sexDim+PA+sol,data=com.dat,lambda="ML");m11.2p<-summary(m11.2)$coefficients["sol",4]
m11.3<-pgls(PA~sexDim+sol,data=com.dat,lambda="ML");m11.3p<-summary(m11.3)$coefficients["sol",4]
m11.4<-pgls(seaDim~sol+sexDim+PA+migDist,data=com.dat,lambda="ML");m11.4p<-summary(m11.4)$coefficients["migDist",4]
m11.5<-pgls(PA~sol+sexDim+migDist,data=com.dat,lambda="ML");m11.5p<-summary(m11.5)$coefficients["migDist",4]

C11<--2*(log(m11.1p)+log(m11.2p)+log(m11.3p)+log(m11.4p)+log(m11.5p))
C11.pval<-1-pchisq(C11,2*5)
CICc11<-CICc(C11,9,108)




 
#Model 12 (3 CIs)
m12.1<-pgls(PA~migDist,data=com.dat,lambda="ML");m12.1p<-summary(m12.1)$coefficients["migDist",4]
m12.2<-pgls(seaDim~migDist,data=com.dat,lambda="ML");m12.2p<-summary(m12.2)$coefficients["migDist",4]
m12.3<-pgls(sexDim~migDist,data=com.dat,lambda="ML");m12.3p<-summary(m12.3)$coefficients["migDist",4]

C12<--2*(log(m12.1p)+log(m12.2p)+log(m12.3p))
C12.pval<-1-pchisq(C12,2*3)
CICc12<-CICc(C12,7,108)




#Model 13 (4 CIs)





m13.1<-pgls(seaDim~sol+PA+migDist,data=com.dat,lambda="ML");m13.1p<-summary(m13.1)$coefficients["migDist",4] 
m13.2<-pgls(sexDim~sol+seaDim+migDist,data=com.dat,lambda="ML");m13.2p<-summary(m13.2)$coefficients["migDist",4]
m13.3<-pgls(sexDim~seaDim+sol,data=com.dat,lambda="ML");m13.3p<-summary(m13.3)$coefficients["sol",4]
m13.4<-pgls(sexDim~sol+seaDim+PA,data=com.dat,lambda="ML");m13.4p<-summary(m13.4)$coefficients["sol",4]

C13<--2*(log(m13.1p)+log(m13.2p)+log(m13.3p)+log(m13.4p))
C13.pval<-1-pchisq(C13,2*4)
CICc13<-CICc(C13,11,108)



#Model 14 (5 CIs)



m14.1<-pgls(PA~sol+migDist,data=com.dat,lambda="ML");m14.1p<-summary(m14.1)$coefficients["migDist",4]
m14.2<-pgls(seaDim~sol+PA+migDist,data=com.dat,lambda="ML");m14.2p<-summary(m14.2)$coefficients["migDist",4]
m14.3<-pgls(migDist~sexDim+sol,data=com.dat,lambda="ML");m14.3p<-summary(m14.3)$coefficients["sexDim",4]
m14.4<-pgls(sexDim~sol,data=com.dat,lambda="ML");m14.4p<-summary(m14.4)$coefficients["sol",4]
m14.5<-pgls(sexDim~sol+PA,data=com.dat,lambda="ML");m14.5p<-summary(m14.5)$coefficients["PA",4]


C14<--2*(log(m1.1p)+log(m1.4p)+log(m14.3p)+log(m14.4p)+log(m14.5p))
C14.pval<-1-pchisq(C14,2*5)
CICc14<-CICc(C14,10,108)


#Model 15




#Model 15 (6 CIs)

#*(mig,pa){dl} 		*pa~dl+mig
#*(mig,td){pa}		*td~pa+mig
#*(dl,td){pa}		td~pa+dl
#(dl,sd){mig}		sd~mig+dl
#(pa,sd){mig,dl}	sd~mig+dl+pa
#(td,sd){pa,mig}	sd~pa+mig+td


m15.1<-m1.1
m15.2<-m1.2
m15.3<-pgls(seaDim~PA+sol,data=com.dat,lambda="ML");m15.3p<-summary(m15.3)$coefficients["sol",4]
m15.4<-pgls(sexDim~migDist+sol,data=com.dat,lambda="ML");m15.4p<-summary(m15.4)$coefficients["sol",4]
m15.5<-pgls(sexDim~migDist+sol+PA,data=com.dat,lambda="ML");m15.5p<-summary(m15.5)$coefficients["PA",4]		
m15.6<-pgls(sexDim~PA+migDist+seaDim,data=com.dat,lambda="ML");m15.6p<-summary(m15.6)$coefficients["seaDim",4]

C15<--2*(log(m1.1p)+log(m1.2p)+log(m1.4p)+log(m15.4p)+log(m15.5p)+log(m15.6p))
C15.pval<-1-pchisq(C15,2*6)
CICc15<-CICc(C15,9,108)



model.all <- c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", 
    "Model 7", "Model 8", "Model 9","Model 10","Model 11","Model 12","Model 13","Model 14","Model 15")
    
C.all <- c(C1, C2, C3, C4, C5, C6, C7, C8, C9,C10,C11,C12,C13,C14,C15)
C.pval.all <- c(C1.pval, C2.pval, C3.pval, C4.pval, C5.pval, C6.pval, C7.pval, C8.pval, C9.pval,C10.pval,C11.pval,C12.pval,C13.pval,C14.pval,C15.pval)
CICc.all <- c(CICc1, CICc2, CICc3, CICc4, CICc5, CICc6, CICc7, CICc8, CICc9, CICc10, CICc11, CICc12, CICc13, CICc14,CICc15)
C.all <- round(C.all, digits = 3)
C.pval.all <- signif(C.pval.all, digits = 3)
CICc.all <- round(CICc.all, digits = 3)
results.dat <- data.frame(model.all, C.all, C.pval.all, CICc.all)
names(results.dat) <- c("Model", "C statistic", "p-value", "CICc")
results.dat <- results.dat[order(results.dat$CICc), ]
print(results.dat) 

write.table(results.dat,file="~/Dropbox/Warbler.Molt.Migration/PATHResults.txt",quote=FALSE)





library(phytools)
library(ggplot2)
library(geiger)

ppTrees<-read.nexus("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Emberizoid_trees/sptree_posterior_4000.nex") 

war_PA<-read.csv("/Users/ryanterrill/Dropbox/Warbler.Molt.Migration/Evolution\ figures/Parulidae_PA.csv")
war_PA[is.na(war_PA)]<-0


notWarbs<-setdiff(ppTrees[[1]]$tip.label,war_PA$tree_taxon)


ppWarb<-llply(ppTrees,drop.tip,tip=notWarbs,.progress="text")
class(ppWarb)<-"multiPhylo"

###Calculate lambda, AIC, logLik, and brownian parameters under ER and ARD models

#ER

pa_ER<-llply(ppWarb,fitDicsrete,
sex_dim_ER
season_dim_ER
migratio_ER
daylength_ER

#ARD

pa_ARD
sex_dim_ARD
season_dim_ARD
migratio_ARD
daylength_ARD
