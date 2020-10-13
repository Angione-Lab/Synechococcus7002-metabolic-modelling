#Set working directory and path
path="C:/Users/Supreeta/Dropbox/final_thesis_Supreeta/"
setwd(path)
getwd()

#Load libraries
library(fmsb)

# Construct the data set
# ATP data
ATP_radar_data <- data.frame(Amino_acid_biosynthesis=c(0.6,0,0.30651242,0.01554119),
Amino_acid_metabolism=c(0.6,0,0.177967687,0.163722174),
Aminoacyl_tRNA_biosynthesis=c(0.6,0,0.30651242,0.015541193),
Biomass_synthesis=c(0.6,0,0.076628105,0.003885299),
Biosynthesis_of_other_secondary_metabolites=c(0.6,0,0.07666572,0.231474948),
#Biotin_biosynthesis=c(0.6,0,0,0),
Carbohydrate_metabolism=c(0.6,0,0.082585496,0.136718442),
Carotenoid_biosynthesis=c(0.6,0,0.203839767,0.010180685),
Cell_wall=c(0.6,0,0.204186696,0.010213233),
Coenzymes_and_prosthetic_groups=c(0.6,0,0.210660378,0.463647872),
#Cyanophycin_metabolism=c(0.6,0,0,0),
Energy_metabolism=c(0.6,0,0.1213941,0.089979245),
Exchange=c(0.6,0,0.114476013,0.242588854),
Fatty_acid_synthesis=c(0.6,0,0.048150944,0.105976656),
Folate_metabolism=c(0.6,0,0.289239578,0.004217668),
Glycan_biosynthesis_and_metabolism=c(0.6,0,0.174917548,0.008660275),
#Glycerophospholipid_metabolism=c(0.6,0,0,0),
Hydrogen_metabolism=c(0.6,0,0.000802431,0.516850126),
Lipid_metabolism=c(0.6,0,0.230057359,0.094139401),
Membrane_bioenergetics=c(0.6,0,0.105330189,0.231823936),
Metabolism_of_cofactors_and_vitamins=c(0.6,0,0.160385831,0.147360916),
Metabolism_of_other_amino_acids=c(0.6,0,0.159884502,0.091611581),
Metabolism_of_terpenoids_and_polyketides=c(0.6,0,0.218863582,0.011068623),
Modelling=c(0.6,0,0.251120231,0.277431343),
None=c(0.6,0,0.025041117,0.090878451),
Nucleotide_metabolism=c(0.6,0,0.133203214,0.157787149),
#Nucleotides_and_nucleic_acids=c(0.6,0,0,0),
Peptidoglycan_biosynthesis=c(0.6,0,0.306047607,0.015100378),
#Phenylmercury_acetate_degradation=c(0.6,0,0,0),
Proline_biosynthesis=c(0.6,0,0.272841115,0.011248198),
#Purine_metabolism=c(0.6,0,0,0),
#Pyridine_metabolism=c(0.6,0,0,0),
#Quinone_biosynthesis=c(0.6,0,0,0),
Riboflavin_metabolism=c(0.6,0,0.236652698,0.353545794),
#Salt_tolerance=c(0.6,0,0,0),
Thiamine_metabolism=c(0.6,0,0.210660378,0.463647872),
Transport=c(0.6,0,0.105736385,0.199163875),
Extracellular_transport=c(0.6,0,0.005438264,0.341758109),
#Vitamin_E_biosynthesis=c(0.6,0,0,0),
row.names = c("max", "min", "Component_1_Average", "Component_2_Average"))

#P1 data
P1_radar_data <- data.frame(Amino_acid_biosynthesis=c(0.6,0,0.26911335,0.006057474),
Amino_acid_metabolism=c(0.6,0,0.16474525,0.049523445),
Aminoacyl_tRNA_biosynthesis=c(0.6,0,0.269184329,0.004497699),
Biomass_synthesis=c(0.6,0,0.067278337,0.001514369),
Biosynthesis_of_other_secondary_metabolites=c(0.6,0,0.130438376,0.016839057),
#Biotin_biosynthesis=c(0.6,0,0,0),
Carbohydrate_metabolism=c(0.6,0,0.101772789,0.109891169),
Carotenoid_biosynthesis=c(0.6,0,0.178951718,0.003869163),
Cell_wall=c(0.6,0,0.179250496,0.002931104),
Coenzymes_and_prosthetic_groups=c(0.6,0,0.179143371,0.378323528),
#Cyanophycin_metabolism=c(0.6,0,0,0),
Energy_metabolism=c(0.6,0,0.157176633,0.226793687),
Exchange=c(0.6,0,0.157877241,0.153980055),
Fatty_acid_synthesis=c(0.6,0,0.037083115,0.111438633),
Folate_metabolism=c(0.6,0,0.25230296,0.0000605305),
Glycan_biosynthesis_and_metabolism=c(0.6,0,0.153567213,0.003387666),
#Glycerophospholipid_metabolism=c(0.6,0,0,0),
Hydrogen_metabolism=c(0.6,0,0.026146192,0.30522191),
Lipid_metabolism=c(0.6,0,0.199762321,0.130447813),
Membrane_bioenergetics=c(0.6,0,0.089571685,0.189161764),
Metabolism_of_cofactors_and_vitamins=c(0.6,0,0.142859159,0.106570622),
Metabolism_of_other_amino_acids=c(0.6,0,0.142536162,0.152387083),
Metabolism_of_terpenoids_and_polyketides=c(0.6,0,0.192155717,0.004303269),
Modelling=c(0.6,0,0.216485445,0.234133319),
None=c(0.6,0,0.054077146,0.027515325),
Nucleotide_metabolism=c(0.6,0,0.111279914,0.312269376),
Nucleotides_and_nucleic_acids=c(0.6,0,0,0),
Peptidoglycan_biosynthesis=c(0.6,0,0.268689662,0.005909977),
#Phenylmercury_acetate_degradation=c(0.6,0,0,0),
Proline_biosynthesis=c(0.6,0,0.0000276151,0.175563802),
#Purine_metabolism=c(0.6,0,0,0),
Pyridine_metabolism=c(0.6,0,0.010870439,0.481184263),
#Quinone_biosynthesis=c(0.6,0,0,0),
Riboflavin_metabolism=c(0.6,0,0.201623297,0.318166462),
#Salt_tolerance=c(0.6,0,0,0),
Thiamine_metabolism=c(0.6,0,0.179143371,0.378323528),
Transport=c(0.6,0,0.139762089,0.25393944),
Extracellular_Transport=c(0.6,0,0.087636239,0.070832862),
#Vitamin_E_biosynthesis=c(0.6,0,0,0),
row.names = c("max", "min", "Component_1_Average", "Component_2_Average"))

#P2 data

P2_radar_data <- data.frame(Amino_acid_biosynthesis=c(1.4,0,0.271600562,0.001347058),
Amino_acid_metabolism=c(1.4,0,0.165695919,0.03889242),
Aminoacyl_tRNA_biosynthesis=c(1.4,0,0.271601482,0.00129183),
Biomass_synthesis=c(1.4,0,0.06790014,0.000336765),
Biosynthesis_of_other_secondary_metabolites=c(1.4,0,0.131830472,0.006718014),
#Biotin_biosynthesis=c(1.4,0,0,0),
Carbohydrate_metabolism=c(1.4,0,0.103150501,0.123538344),
Carotenoid_biosynthesis=c(1.4,0,0.180579434,0.000687324),
Cell_wall=c(1.4,0,0.162016006,0.096166753),
Coenzymes_and_prosthetic_groups=c(1.4,0,0.177724648,0.267076336),
#Cyanophycin_metabolism=c(1.4,0,0,0),
Energy_metabolism=c(1.4,0,0.156674415,0.08766789),
Exchange=c(1.4,0,0.157955965,0.070165372),
Fatty_acid_synthesis=c(1.4,0,0.03640205,0.094706674),
Folate_metabolism=c(1.4,0,0.254076949,0.001405695),
Glycan_biosynthesis_and_metabolism=c(1.4,0,0.154974302,0.000666807),
#Glycerophospholipid_metabolism=c(1.4,0,0,0),
Hydrogen_metabolism=c(1.4,0,0.026892507,0.135827961),
Lipid_metabolism=c(1.4,0,0.203987871,0.115914435),
Membrane_bioenergetics=c(1.4,0,0.088862324,0.133538168),
Metabolism_of_cofactors_and_vitamins=c(1.4,0,0.142968857,0.087654732),
Metabolism_of_other_amino_acids=c(1.4,0,0.137867862,0.068482124),
Metabolism_of_terpenoids_and_polyketides=c(1.4,0,0.193927754,0.000929095),
Modelling=c(1.4,0,0.216490756,0.166754933),
None=c(1.4,0,0.071278835,0.395591186),
Nucleotide_metabolism=c(1.4,0,0.105853617,0.403974347),
#Nucleotides_and_nucleic_acids=c(1.4,0,0,0),
Peptidoglycan_biosynthesis=c(1.4,0,0.271148524,0.001141177),
#Phenylmercury_acetate_degradation=c(1.4,0,0,0),
Proline_biosynthesis=c(1.4,0,0.0000552831,0.03263899),
Purine_metabolism=c(1.4,0,0.03959598,1.361640453),
#Pyridine_metabolism=c(1.4,0,0,0),
#Quinone_biosynthesis=c(1.4,0,0,0),
Riboflavin_metabolism=c(1.4,0,0.200599091,0.232047287),
#Salt_tolerance=c(1.4,0,0,0),
Thiamine_metabolism=c(1.4,0,0.177724648,0.267076336),
Transport=c(1.4,0,0.133608278,0.33261616),
Extracellular_Transport=c(1.4,0,0.088319861,0.044650697),
#Vitamin_E_biosynthesis=c(1.4,0,0,0),
row.names = c("max", "min", "Component_1_Average", "Component_2_Average"))

#Gene transcripts
Gene_radar_data <- data.frame(Amino_acid_transport_and_metabolism=c(0.06,0,0.042236232,0.03372008),
Carbohydrate_transport_and_metabolism=c(0.06,0,0.035161074,0.029810811),
Cell_cycle_control_cell_division_and_chromosome_partitioning=c(0.06,0,0.036092629,0.031758665),
Cell_wall_membrane_and_envelope_biogenesis=c(0.06,0,0.044596947,0.02478318),
Cell_motility=c(0.06,0,0.040175848,0.033622006),
Chromatin_structure_and_dynamics=c(0.06,0,0.019382275,0.059062088),
Coenzyme_transport_and_metabolism=c(0.06,0,0.037768258,0.037559198),
Defense_mechanisms=c(0.06,0,0.030115869,0.027680755),
Energy_production_and_conversion=c(0.06,0,0.033350883,0.033654175),
General_function_prediction_only=c(0.06,0,0.037159566,0.02970199),
Inorganic_ion_transport_and_metabolism=c(0.06,0,0.026378391,0.025252386),
Intracellular_trafficking_secretion_and_vesicular_transport=c(0.06,0,0.035768277,0.039503638),
Lipid_transport_and_metabolism=c(0.06,0,0.035451112,0.036481235),
None=c(0.06,0,0.024976774,0.030379754),
Nucleotide_transport_and_metabolism=c(0.06,0,0.040799353,0.025993845),
Posttranslational_modification_protein_turnover_chaperones=c(0.06,0,0.035657846,0.023061427),
Replication_recombination_and_repair=c(0.06,0,0.038118448,0.027117483),
Secondary_metabolites_biosynthesis_transport_and_metabolism=c(0.06,0,0.028606031,0.031231554),
Signal_transduction_mechanisms=c(0.06,0,0.034901441,0.033039929),
Transcription=c(0.06,0,0.0347851,0.028393347),
Translation_ribosomal_structure_and_biogenesis=c(0.06,0,0.032093636,0.038858733),
row.names = c("max", "min", "Component_1_Average", "Component_2_Average"))

# Define line colors
colors_line_ATP <- c(scales::alpha("green3", 0.9),
                scales::alpha("orangered", 0.9))
colors_line_P1 <- c(scales::alpha("darkblue", 0.9),
                scales::alpha("red3", 0.9))
colors_line_P2 <- c(scales::alpha("gold", 0.9),
                scales::alpha("orchid4", 0.9))
colors_line_gene <- c(scales::alpha("cyan4", 0.9),
                scales::alpha("indianred1", 0.9))

# Create plot
radarchart(ATP_radar_data, 
           seg = 6,  # Number of axis segments
           title = "Average Component Score Radar Chart",
           pcol = colors_line_ATP,
           plty = 1:1,
           plwd = 2,
	     axistype = 4,
           caxislabels = c("0","0.1","0.2","0.3","0.4","0.5","0.6"),
           cglty = 3,
	     cglcol = "gray70",
	     axislabcol="gray0")
# Add a legend
legend(x=1.35, y=1.25, legend = rownames(ATP_radar_data[-c(1,2),]), bty = "o", pch=20 , col=colors_line_ATP , text.col = "gray0", cex=1.2, pt.cex=3)

radarchart(P1_radar_data, 
           seg = 6,  # Number of axis segments
           title = "Average Component Score Radar Chart",
           pcol = colors_line_P1,
	     plty = 1:1,
           plwd = 2,
	     axistype = 4,
	     caxislabels = c("0","0.1","0.2","0.3","0.4","0.5","0.6"),
           cglty = 3,
	     cglcol = "gray70",
	     axislabcol="gray0")
# Add a legend
legend(x=1.35, y=1.25, legend = rownames(P1_radar_data[-c(1,2),]), bty = "o", pch=20 , col=colors_line_P1 , text.col = "gray0", cex=1.2, pt.cex=3)

radarchart(P2_radar_data, 
           seg = 7,  # Number of axis segments
           title = "Average Component Score Radar Chart",
           pcol = colors_line_P2,
	     plty = 1:1,
           plwd = 2,
	     axistype = 4,
	     caxislabels = c("0","0.2","0.4","0.6","0.8","1.0","1.2","1.4"),
           cglty = 3,
	     cglcol = "gray70",
	     axislabcol="gray0")
# Add a legend
legend(x=1.35, y=1.25, legend = rownames(P2_radar_data[-c(1,2),]), bty = "o", pch=20 , col=colors_line_P2 , text.col = "gray0", cex=1.2, pt.cex=3)

radarchart(Gene_radar_data, 
           seg = 6,  # Number of axis segments
           title = "Average Component Score Radar Chart",
           pcol = colors_line_gene,
	     plty = 1:1,
           plwd = 2,
	     axistype = 4,
	     caxislabels = c("0","0.01","0.02","0.03","0.04","0.05","0.06"),
           cglty = 3,
	     cglcol = "gray70",
	     axislabcol="gray0")
# Add a legend
legend(x=1.35, y=1.25, legend = rownames(Gene_radar_data[-c(1,2),]), bty = "o", pch=20 , col=colors_line_gene , text.col = "gray0", cex=1.2, pt.cex=3)





