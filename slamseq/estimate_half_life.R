library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

give.n <- function(x){
    return(c(y = median(x)*1.05, label = length(x)))
}

corticali3_gran = ReadGRAND('exp1_grandslam_grandR.tsv',
                      design=c(Design$Condition,
                               Design$Replicate,
                               Design$dur.4sU))


new = Coldata(corticali3_gran) |> 
    mutate(timeCondition = paste0(Condition,duration.4sU)) |> 
    mutate(timeCondition = as.factor(timeCondition)) |> 
    mutate(no4sU = ifelse(duration.4sU.original == "0h", TRUE,FALSE ))

corticali3_gran$coldata = new


corticali3_gran = FilterGenes(corticali3_gran)
PlotPCA(corticali3_gran,aest=aes(color=Condition,shape=as.factor(duration.4sU)))

# kinetic modelling  -------------------------------------------------

corticali3_norm <-Normalize(corticali3_gran)


SetParallel(cores = 4)  # increase 2 on your system, or omit the cores = 2 for automatic detection

corticali3_norm<-FitKinetics(corticali3_norm,name = "kinetics")

# exp2_norm<-FitKinetics(exp2_norm,name = "kinetics")
# exp2_kinetic <-GetAnalysisTable(exp2_norm)

corticali3_kinetic <-GetAnalysisTable(corticali3_norm)

corticali3_kinetic = corticali3_kinetic|> 
    as.data.frame() |> 
    tibble::rownames_to_column('gene_name') |> 
    as.data.table()

# exp2_kinetic = exp2_kinetic |> as.data.frame() |> 
#     tibble::rownames_to_column('gene_name') |> 
#     as.data.table()

corticali3_kinetic = corticali3_kinetic |> 
    mutate(log2FoldHalfLife = log2(`kinetics.TDPKD.Half-life` /`kinetics.Ctrl.Half-life` ))

# exp2_kinetic = exp2_kinetic |> mutate(log2FoldHalfLife = log2(`kinetics.TDP43kd.Half-life` /`kinetics.control.Half-life` ))

PlotScatter(corticali3_kinetic,
            analysis="kinetics",
            x=`kinetics.Ctrl.Half-life`,
            y=`kinetics.TDPKD.Half-life`,log=TRUE)+
    geom_abline()



corticali3_kinetic |> 
    select(`kinetics.Ctrl.Half-life`,`kinetics.TDPKD.Half-life`,Gene) |> 
    melt() |> filter(value < 30) |> 
    ggplot(aes(y = value, x = variable)) + 
    geom_boxplot() +
    ggpubr::stat_compare_means() +
    stat_summary(fun.data = give.n, geom = "text", fun.y = median,
                 position = position_dodge(width = 0.75)) + 
    ylim(0,24)


exp1_kinetic |> 
    filter(`kinetics.Ctrl.Half-life` < 40 & `kinetics.Ctrl.Half-life` > 2) |> 
    mutate(binned_control_half = cut_number(`kinetics.Ctrl.Half-life`,n = 5)) |> 
    ggplot(aes(x =binned_control_half,y =  log2FoldHalfLife)) + 
    geom_boxplot()



# Generate the half-life plots ----------------------------------------------

exp1_norm = ComputeNtrCI(exp1_norm)

all_old_rna_long_real = data.table()
all_old_rna_long_fitted = data.table()

this_bois = c("ELK1","SIX3","TLX1")    
for(g in this_bois){
    
    e = PlotGeneProgressiveTimecourse(exp1_norm,g,show.CI = TRUE,return.tables = TRUE)
    old_rna_fit = e$fitted |> 
        filter(Type == 'Old') |> 
        mutate(gene = g)
    all_old_rna_long_fitted = rbind(all_old_rna_long_fitted,old_rna_fit)
    
    old_real = e$df |> 
        filter(Type == 'Old') |> 
        mutate(gene = g)
    
    all_old_rna_long_real = rbind(all_old_rna_long_real, old_real)
}



go_of_interest = 'ELK1'

half_life_control_estimate  = exp1_kinetic |> filter(gene_name == go_of_interest) |> 
    pull(`kinetics.Ctrl.Half-life`)

half_life_tdpkd_estimate = exp1_kinetic |> filter(gene_name == go_of_interest) |> 
    pull(`kinetics.TDPKD.Half-life`)

ggplot() + 
    geom_line(aes(x = time,y = Value,color = Condition),data = all_old_rna_long_fitted,size = 2) + 
    geom_point(aes(x = time,y = Value,color = Condition),data = all_old_rna_long_real,show_guide = FALSE) + 
    geom_pointrange(aes(x = time,color = Condition,ymin=lower, ymax=upper,y = Value),data = all_old_rna_long_real,show_guide = FALSE) + 
    ggpubr::theme_pubr() + 
    ylab("Unlabelled RNA") + 
    xlab("4SU labelling time") + 
    scale_color_manual(values=c('#999999','#E69F00')) + 
    annotate("text",  x=10, y = 1000, label = glue::glue("Control Half-Life  {round(half_life_control_estimate,1)} h"),vjust=1, hjust=0) +
    annotate("text",  x=10, y = 900, label = glue::glue("TDP-43 KD Half-Life {round(half_life_tdpkd_estimate,1)} h"),vjust=1, hjust=0) +
    ggpubr::theme_pubr() +
    ggpubr::labs_pubr() +
    theme(legend.title = element_blank()) +
    scale_x_continuous(breaks = c(0, 1, 4, 8, 12, 24)) + 
    facet_wrap(~gene)
