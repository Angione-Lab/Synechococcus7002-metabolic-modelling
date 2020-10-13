#Pyramid charts of component sum averages by pathways/COG categories
#Set working directory
path="C:/Users/Supreeta/Dropbox/final_thesis_Supreeta/"
setwd(path)
getwd()

#Plotrix's pyramid.plot
#Load library
library(plotrix)

#(Each variable is separate and read in reverse so categories are plotted the right way up)
#Load pathway labels 
pathways<-c("Vitamin E biosynthesis","Extracellular transport","Transport","Thiamine metabolism","Salt tolerance",
"Riboflavin metabolism","Quinone biosynthesis","Pyridine metabolism","Purine metabolism","Proline biosynthesis",
"Phenylmercury acetate degradation","Peptidoglycan biosynthesis","Nucleotides and nucleic acids","Nucleotide metabolism",
"Unassigned","Modelling","Metabolism of terpenoids and polyketides","Metabolism of other amino acids",
"Metabolism of cofactors and vitamins","Membrane bioenergetics","Lipid metabolism","Hydrogen metabolism",
"Glycerophospholipid metabolism","Glycan biosynthesis and metabolism","Folate metabolism","Fatty acid synthesis",
"Exchange","Energy metabolism","Cyanophycin metabolism","Coenzymes and prosthetic groups","Cell wall","Carotenoid biosynthesis",
"Carbohydrate metabolism","Biotin biosynthesis","Biosynthesis of other secondary metabolites","Biomass synthesis",
"Aminoacyl-tRNA biosynthesis","Amino acid metabolism","Amino acid biosynthesis")

#Load ATP data for Component 1
comp1atp.pop<-c(0,0.065259164,5.392555628,0.210660378,0,0.236652698,0,0,0,0.272841115,0,0.306047607,0,10.65625716,
0.450740112,0.753360692,4.596135221,3.037805537,18.44437061,0.210660378,9.432351729,0.001604862,0,2.448845667,
0.289239578,1.685283024,6.296180727,7.769222425,0,0.210660378,0.612560087,4.892154412,7.680451156,0,0.306662878,
0.30651242,0.61302484,24.91547611,0.30651242)
#Load ATP data for Component 2
comp2atp.pop<-c(0,4.101097313,10.15735761,0.463647872,0,0.353545794,0,0,0,0.011248198,0,0.015100378,0,12.6229719,
1.635812122,0.832294029,0.232441083,1.740620034,16.94650538,0.463647872,3.85971545,1.033700252,0,0.121243849,
0.004217668,3.709182976,13.34238699,5.758671693,0,0.463647872,0.0306397,0.244336441,12.71481512,0,0.92589979,
0.015541195,0.031082385,22.92110436,0.01554119)

#Set ATP color gradient using preset-color-palettes from R-colorspace
library("colorspace")
comp1atpcol<-sequential_hcl(9, "Greens")
comp2atpcol<-sequential_hcl(9, "Oranges")

#Plot ATP pyramid
par(mar=pyramid.plot(comp1atp.pop,comp2atp.pop,labels=pathways,
main="Biomass-ATP Component Sum",top.labels=c("Component 1 Sum","Pathway","Component 2 Sum"),unit="", lxcol=comp1atpcol,rxcol=comp2atpcol,
gap=0,xlim=c(25,25),show.values=FALSE))

#Load P1 data for Component 1
comp1p1.pop<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,1.0275228,0,0,0,1.57975936,2.449298591,0.065222632,1.970077836,0.173411093,
0.094159381,2.196799519,0.16931122,6.295673624,11.95276816,16.13981602,1.034491933,0.259329296,0.392806369,6.362411294,
24.96149716,0.269107083,0.807333781,1.076452197,0.269113054,36.77626825,0.270059924)
#Load P1 data for Component 2
comp2p1.pop<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0.026440752,0,0,0,0.038579476,31.96088506,2.88710558,48.40913893,1.193228115,
0.571154617,1.159587452,0.000549402,11.9521439,11.15380476,0.61832349,0.00498442,0.010457692,0.302122598,0.193740339,
0.540472756,0.006053828,0.018168789,0.02422931,0.006057334,0.776046384,0.000534213)

#Set P1 color gradient using preset-color-palettes from  R-colorspace
comp1p1col<-sequential_hcl(9, "Teal")
comp2p1col<-sequential_hcl(9, "Reds")

#Plot P1 pyramid
par(mar=pyramid.plot(comp1p1.pop,comp2p1.pop,labels=pathways,
main="Biomass-Photosystem 1 Component Sum",top.labels=c("Component 1 Sum","Pathway","Component 2 Sum"),unit="", lxcol=comp1p1col,rxcol=comp2p1col,
gap=0,xlim=c(50,50),show.values=FALSE))

#Load P2 data for Component 1
comp1p2.pop<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,1.035404336,0,0,0,1.592787681,2.592283919,0.079191486,1.934141502,
0.135178416,0.067589208,1.94995157,0.140219386,6.059221767,11.83986976,16.24057951,1.042000839,0.260896089,
0.358364854,6.40857835,25.18811477,0.271431092,0.814626519,1.086251899,0.271581248,37.10905196,0.272314341)

#Load P2 data for Component 2
comp2p2.pop<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0.020366378,0,0,0,0.019565956,40.62872933,2.719960333,47.21304775,
0.105979251,0.052989625,2.704680932,0.03852745,7.939708713,8.579536971,0.493892082,0.006547474,0.002280216,
1.564416698,0.154910154,0.123340693,0.001269793,0.003270225,0.001076603,0.000449333,0.202916674,0.00001359)

#Set P2 color gradient using preset-color-palettes from  R-colorspace
comp1p2col<-sequential_hcl(9, "YlOrRd")
comp2p2col<-sequential_hcl(9, "PurpOr")

#Plot P2 pyramid
par(mar=pyramid.plot(comp1p2.pop,comp2p2.pop,labels=pathways,
main="Biomass-Photosystem 2 Component Sum",top.labels=c("Component 1 Sum","Pathway","Component 2 Sum"),unit="",lxcol=comp1p2col,rxcol=comp2p2col,
gap=0,xlim=c(50,50),show.values=FALSE))

#Load COG category labels in reverse
COGpathways<-c("Translation, ribosomal structure and biogenesis","Transcription","Signal transduction mechanisms",
"Secondary metabolites biosynthesis, transport and metabolism","Replication, recombination, and repair",
"Posttranslational modification, protein turnover, chaperones","Nucleotide transport and metabolism",
"Unassigned","Lipid transport and metabolism","Intracellular trafficking, secretion, and vesicular transport",
"Inorganic ion transport and metabolism","General","Energy production and conversion","Defense mechanisms",
"Coenzyme transport and metabolism","Chromatin structure and dynamics","Cell wall, membrane, envelope biogenesis",
"Cell motility","Cell cycle control, cell division, chromosome partitioning",
"Carbohydrate transport and metabolism","Amino acid transport and metabolism")

##Load transcript data for Component 1
comp1gene.pop<-c(8.82574991,4.03507162,5.72383640,1.34448345,4.49797688,4.27894153,
2.44796117,30.87129278,1.84345781,2.07456008,4.37881293,11.07355053,4.86922895,
1.59614108,5.09871478,0.03876455,6.73413900,1.20527545,1.11887149,3.30514097,7.09568691)

#Load transcript data for Component 2
comp2gene.pop<-c(10.68615163,3.29362826,5.41854833,1.46788302,3.19986297,2.76737124,
1.55963069,37.54937544,1.89702420,2.29121103,4.19189607,8.85119309,4.91350962,
1.46708003,5.07049174,0.11812418,3.74226012,1.00866017,0.98451863,2.80221622,5.66497351)

#Set transcript color gradient using preset-color-palettes from  R-colorspace
comp1genecol<-sequential_hcl(9, "BluGrn")
comp2genecol<-sequential_hcl(9, "Peach")

#Plot transcript pyramid
par(mar=pyramid.plot(comp1gene.pop,comp2gene.pop,labels=COGpathways,
main="Gene Transcript Component Sum",top.labels=c("Component 1 Sum","COG Category","Component 2 Sum"),unit="",lxcol=comp1genecol,rxcol=comp2genecol,
gap=0,xlim=c(40,40),show.values=FALSE))