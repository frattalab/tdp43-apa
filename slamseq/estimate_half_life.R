library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(grandR)
library(ggrastr)
library(glue)
give.n <- function(x){
    return(c(y = median(x)*1.05, label = length(x)))
}

#' function to plot fitted half life curves in two conditions for specific genes
plot_hl <- function(go_of_interest, fitted = all_old_rna_long_fitted, real = all_old_rna_long_real, kinetic = corticali3_kinetic) {
  
  # extract CTL & KD half lives for gene of interest
  go_kinetic <- filter(kinetic, Symbol == go_of_interest)
  
  ctl_hl <- pull(go_kinetic, `kinetics.Ctrl.Half-life`)
  kd_hl <- pull(go_kinetic, `kinetics.TDPKD.Half-life`)
  
  ggplot() + 
    geom_line(aes(x = time,y = Value,color = Condition),data = all_old_rna_long_fitted[gene == go_of_interest],linewidth = 1.5) + 
    geom_point(aes(x = time,y = Value,color = Condition),data = all_old_rna_long_real[gene == go_of_interest], show.legend =  FALSE) + 
    geom_pointrange(aes(x = time,color = Condition,ymin=lower, ymax=upper,y = Value),data = all_old_rna_long_real[gene == go_of_interest],
                    show.legend = FALSE) + 
    ylab("Unlabelled RNA") + 
    xlab("4SU labelling time") + 
    scale_color_manual(values=c('#999999','#E69F00')) +
    scale_x_continuous(breaks = c(0, 1, 4, 8, 12, 24)) + 
    scale_y_continuous(limits = c(0,1500),
                       breaks = seq(0,1500,250)) +
    geom_text(aes(x = 24, y = 1500, label = glue("Control: {round(ctl_hl,2)} hr")),
              hjust = 1, #right align
              size = 7.5, colour = "#999999") +
    geom_text(aes(x = 24, y = 1400, label = glue("TDP-43 KD: {round(kd_hl,2)} hr")),
              hjust = 1, #right align
              size = 7.5, colour = "#E69F00") +
    theme_classic(base_size = 20) +
    # ggpubr::theme_pubr(base_size = 20) +
    theme(legend.title = element_blank(),
          legend.position = "top") +
    facet_wrap(~gene)
  
}


corticali3_gran = ReadGRAND('data/grandslam_exp1_grandR.tsv',
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


corticali3_kinetic = corticali3_kinetic |> 
    mutate(log2FoldHalfLife = log2(`kinetics.TDPKD.Half-life` /`kinetics.Ctrl.Half-life` ))

# Also calculate change in synthesis rate
corticali3_kinetic = corticali3_kinetic |>
  mutate(log2FoldSynthesis = log2(kinetics.TDPKD.Synthesis / kinetics.Ctrl.Synthesis))

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


corticali3_kinetic |> 
    filter(`kinetics.Ctrl.Half-life` < 40 & `kinetics.Ctrl.Half-life` > 2) |> 
    mutate(binned_control_half = cut_number(`kinetics.Ctrl.Half-life`,n = 5)) |> 
    ggplot(aes(x =binned_control_half,y =  log2FoldHalfLife)) + 
    geom_boxplot()



# Generate the half-life plots ----------------------------------------------

corticali3_norm = ComputeNtrCI(corticali3_norm)

all_old_rna_long_real = data.table()
all_old_rna_long_fitted = data.table()

this_bois = c("ELK1","SIX3","TLX1", "BRINP2")    

for(g in this_bois){
    
    e = PlotGeneProgressiveTimecourse(corticali3_norm,g,show.CI = TRUE,return.tables = TRUE)
    old_rna_fit = e$fitted |> 
        filter(Type == 'Old') |> 
        mutate(gene = g)
    all_old_rna_long_fitted = rbind(all_old_rna_long_fitted,old_rna_fit)
    
    old_real = e$df |> 
        filter(Type == 'Old') |> 
        mutate(gene = g)
    
    all_old_rna_long_real = rbind(all_old_rna_long_real, old_real)
}



genes_of_interest <- this_bois %>% set_names()

# plots for each 3'UTR cryptic
goi_hl_curves <- map(genes_of_interest,
    ~ plot_hl(.x))


if (!dir.exists("processed")) {dir.create("processed", recursive = T)}

fwrite(corticali3_kinetic, "processed/2023-08-22_i3cortical_slamseq_grandr_kinetics.tsv",sep = "\t",col.names = T)

walk2(goi_hl_curves,
      names(goi_hl_curves),
      ~ ggsave(filename = glue("2023-10-10_i3cortical_hl_curve_{.y}.svg"),
               plot = .x,
               device = svg,
               path = "processed/",
               width = 6,
               height = 6,
               units = "in",
               dpi = "retina"
               )
      )
