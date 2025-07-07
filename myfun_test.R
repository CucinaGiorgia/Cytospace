#LOADING RAW DATA

loading_data <- function(file_path){
  raw_data <- data.table::fread(file_path) %>% #Carica file
    tibble::as_tibble(.name_repair=janitor::make_clean_names) #Modifica file: il nome delle colonne senza spazi e con minuscolo
  return(raw_data)
}

#MAKE EXPERIMENTAL DESIGN
exp_design <- function(data, pattern_interest){
  
  exp_des <- data %>% 
    dplyr::select(gene_names, dplyr::starts_with(pattern_interest)) %>%
    tidyr::pivot_longer(!gene_names, names_to="key", values_to= "intensity") %>% #Raggruppa tutte le lfq in unica colonna
    dplyr::distinct(key) %>%
    dplyr::mutate(label = stringr::str_remove(key, "lfq_intensity_")) %>%
    dplyr::mutate(condition = stringr::str_remove(label, "_[^_]*$")) %>%
    dplyr::mutate(replicate = stringr::str_remove(label, ".*_"))
  
  return(exp_des)
}

#DATA PREPROCESSING
pre_process <- function(data, pattern_interest){
  
  data_pre <- data %>% 
    dplyr::mutate(dplyr::across(dplyr::starts_with(pattern_interest), ~ log2(.))) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with(pattern_interest), ~ dplyr::na_if(.,-Inf)))
  
  return(data_pre)
}

#DATA WRANGLING
data_wrangling <- function(data, 
                           pep_filter, 
                           pep_thr, 
                           rev=TRUE,
                           cont=TRUE,
                           oibs=TRUE){
  
  data_wrang <- data %>% #DATA STANDARDIZED
    dplyr::select(protein_i_ds, gene_names, id) %>%
    dplyr::mutate(gene_names = stringr::str_extract(gene_names, "[^;]*")) %>%
    ## every protein groups now have only 1 gene name associated to it
    dplyr::rename(unique_gene_names = gene_names) %>%
    janitor::get_dupes(unique_gene_names) %>%
    dplyr::mutate(unique_gene_names = dplyr::case_when(
      unique_gene_names != "" ~ paste0(
        unique_gene_names, "__",
        stringr::str_extract(protein_i_ds, "[^;]*")),
      TRUE ~ stringr::str_extract(protein_i_ds, "[^;]*"))) %>%
    dplyr::select(unique_gene_names, id) %>%
    dplyr::right_join(data, by = "id") %>%
    dplyr::mutate(gene_names = dplyr::case_when(unique_gene_names != "" ~ unique_gene_names,
                                                TRUE ~ gene_names)) %>%
    dplyr::select(-unique_gene_names) %>%
    dplyr::mutate(gene_names = dplyr::if_else(gene_names == "",
                                              stringr::str_extract(protein_i_ds, "[^;]*"),
                                              gene_names)) %>%
    dplyr::mutate(gene_names = stringr::str_extract(gene_names, "[^;]*")) %>% 
    dplyr::select(gene_names,
                  dplyr::all_of(exp_des$key),
                  peptides,
                  razor_unique_peptides,
                  unique_peptides,
                  reverse,
                  potential_contaminant,
                  only_identified_by_site) %>% 
    tidyr::pivot_longer(!c(gene_names,
                           peptides,
                           razor_unique_peptides,
                           unique_peptides,
                           reverse,
                           potential_contaminant,
                           only_identified_by_site),
                        names_to = "key",
                        values_to = "raw_intensity") %>% 
    dplyr::inner_join(., exp_des, by = "key") %>%   #aggiunge righe e colonne che matchano tra expdesign e data
    dplyr::mutate(bin_intensity = dplyr::if_else(is.na(raw_intensity), 0, 1)) %>%  #Nuova colonna. 1 se valore esite, 0 se NA
    dplyr::select(-key) %>% 
    {if(rev)dplyr::filter(., !reverse == "+") else .} %>% #DATA WRANGLING
    {if(cont)dplyr::filter(., !potential_contaminant == "+") else .} %>%
    {if(oibs)dplyr::filter(., !only_identified_by_site == "+") else .} %>% 
    ## filter on peptides:
    {if(pep_filter == "peptides"){dplyr::filter(., peptides >= pep_thr)}
      else if (pep_filter == "unique") {dplyr::filter(., unique_peptides >= pep_thr)}
      else {dplyr::filter(., razor_unique_peptides >= pep_thr)}}
  
  return(data_wrang)
  
}

#DATA FILTERING
data_filtered <- function(data,
                          valid_val_filter,
                          valid_val_thr) {
  data_filt<- data %>% 
    {if(valid_val_filter == "total")dplyr::group_by(., gene_names)
      else dplyr::group_by(., gene_names, condition)} %>% 
    dplyr::mutate(miss_val = dplyr::n() - sum(bin_intensity)) %>% 
    dplyr::mutate(n_size = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(gene_names) %>%
    ## rage compreso tra 0 e 100% espresso in valori tra 0 e 1
    {if(valid_val_filter == "alog") dplyr::filter(., any(miss_val <= round(n_size * (1 - valid_val_thr), 0)))
      else dplyr::filter(., all(miss_val <= round(n_size * (1 - valid_val_thr), 0)))} %>%
    dplyr::ungroup() %>%
    dplyr::select(gene_names, label, condition, replicate, bin_intensity, raw_intensity) %>% 
    dplyr::rename(intensity = raw_intensity)
  return(data_filt)
}

#DATA IMPUTATION
data_imputed <- function(data,
                         shift, 
                         scale, 
                         unique_visual = FALSE){
  
  imputed_data <- data %>%
    dplyr::group_by(gene_names, condition) %>%
    dplyr::mutate(for_mean_imp = dplyr::if_else((sum(bin_intensity) / dplyr::n()) >= 0.75, TRUE, FALSE)) %>%
    dplyr::mutate(mean_grp = mean(intensity, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(imp_intensity = dplyr::case_when(
      bin_intensity == 0 & for_mean_imp ~ mean_grp,
      TRUE ~ as.numeric(intensity))) %>%
    dplyr::mutate(intensity = imp_intensity) %>% 
    dplyr::select(-c(for_mean_imp, mean_grp, imp_intensity))%>%
    dplyr::group_by(label) %>%
    # Define statistic to generate the random distribution relative to sample
    dplyr::mutate(mean = mean(intensity, na.rm = TRUE),
                  sd = sd(intensity, na.rm = TRUE),
                  n = sum(!is.na(intensity)),
                  total = nrow(data) - n) %>%
    dplyr::ungroup() %>%
    # Impute missing values by random draws from a distribution
    # which is left-shifted by parameter 'shift' * sd and scaled by parameter 'scale' * sd.
    dplyr::mutate(imp_intensity = dplyr::case_when(is.na(intensity) ~ rnorm(total,
                                                                            mean = (mean - shift * sd), 
                                                                            sd = sd * scale),
                                                   TRUE ~ intensity)) %>%
    dplyr::mutate(intensity = imp_intensity) %>%
    dplyr::select(-c(mean, sd, n, total, imp_intensity)) %>%
    dplyr::group_by(condition)
  
  return(imputed_data) 
}

#DEFINISCI CONFRONTI PER TEST STATISTICO
define_tests <- function(){
  conditions <-
    dplyr::distinct(exp_des, condition) %>% pull(condition)
  
  tests <-
    tidyr::expand_grid(cond1 = conditions, cond2 = conditions) %>%
    dplyr::filter(cond1 != cond2) %>%
    dplyr::mutate(test = paste0(cond1, "_vs_", cond2)) %>%
    dplyr::pull(test)
  
  return(tests)
}

#ANALISI STATISTICA T-TEST
stat_t_test_single <- function(data, test, fc, alpha, p_adj_method, paired_test){
  
  cond_1 <- stringr::str_split(test, "_vs_")[[1]][1]
  cond_2 <- stringr::str_split(test, "_vs_")[[1]][2]
  
  mat <- data %>%
    dplyr::filter(condition == cond_1 | condition == cond_2) %>%
    dplyr::mutate(label_test = paste(condition, replicate, sep = "_")) %>%
    tidyr::pivot_wider(id_cols = "gene_names",
                       names_from = "label_test",
                       values_from = "intensity") %>%
    tibble::column_to_rownames("gene_names") %>%
    dplyr::relocate(dplyr::contains(cond_2), .after = dplyr::last_col()) %>%
    na.omit() %>% 
    as.matrix()
  
  a <- grep(cond_1, colnames(mat))
  b <- grep(cond_2, colnames(mat))
  
  p_values_vec <- apply(mat, 1, function(x) t.test(x[a], x[b], paired_test=FALSE, var_equal=TRUE)$p.value)
  
  p_values <- p_values_vec %>%
    as_tibble(rownames = NA) %>%
    tibble::rownames_to_column(var = "gene_names") %>%
    dplyr::rename(p_val = value)
  
  fold_change <- apply(mat, 1, function(x) mean(x[a]) - mean(x[b])) %>% #metterlo in log2?
    as_tibble(rownames = NA) %>%
    tibble::rownames_to_column(var = "gene_names") %>%
    dplyr::rename(fold_change = value)
  
  p_adjusted <- p.adjust(p_values_vec, method = p_adj_method) %>% 
    as_tibble(rownames = NA) %>%
    tibble::rownames_to_column(var = "gene_names") %>%
    dplyr::rename(p_adj = value)
  
  stat_data <- fold_change %>% 
    dplyr::full_join(., p_values, by = "gene_names") %>% 
    dplyr::full_join(., p_adjusted, by = "gene_names") %>% 
    dplyr::mutate(significant = dplyr::if_else(abs(fold_change) >= fc & p_adj <= alpha, TRUE, FALSE)) %>% 
    dplyr::rename(!!paste0(cond_1, "_vs_", cond_2, "_significant") := significant) %>% 
    dplyr::rename(!!paste0(cond_1, "_vs_", cond_2, "_p_val") := p_val) %>% 
    dplyr::rename(!!paste0(cond_1, "_vs_", cond_2, "_fold_change") := fold_change) %>% 
    dplyr::rename(!!paste0(cond_1, "_vs_", cond_2, "_p_adj") := p_adj)
  
  return(stat_data)
}

#COLOR PALETTE (8 options, from A to h)
define_colors = function(){
  n_of_color <- max(exp_des %>% dplyr::count(replicate) %>% dplyr::pull(n))
  color_palette <- viridis::viridis(n = n_of_color , direction = -1, end = 0.70, begin = 0.30)
}

#BARPLOT COUNT
barplot_count <- function (data,
                           color="black",
                           width_barplot=0.5){
  bar <- data %>% 
    dplyr::group_by(label) %>%
    dplyr::summarise(counts = sum(bin_intensity)) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(., exp_des, by = "label") %>%
    dplyr::mutate(replicate = as.factor(replicate)) %>%
    dplyr::group_by(condition) %>% 
    
    ggplot2::ggplot(aes(x=label, y=counts, fill=condition))+
    geom_bar(stat="identity", width = width_barplot, color=color)+
    scale_fill_manual(values=color_palette)+
    theme_cuc()+
    theme(axis.text.x = element_text(angle = 30))+
    labs(title= "Protein counts", x="Label", y="Counts")+
    geom_text(aes(label=counts), vjust=-0.2, size=4)
  
  return(bar)
}

#BARPLOT COVERAGE
barplot_cover <- function (data,
                           fillcolor="#21908CFF",
                           color="black",
                           width_barplot=0.5){
  bar <- data %>% 
    dplyr::group_by(gene_names) %>%
    dplyr::summarise(counts = sum(bin_intensity)) %>%
    dplyr::ungroup() %>%
    dplyr::count(counts) %>% 
    dplyr::rename(occurrence = n) %>% 
    
    ggplot2::ggplot(aes(x=counts, y=occurrence))+
    geom_bar(stat="identity", width = width_barplot, color=color, fill=fillcolor)+
    theme_cuc()+
    scale_fill_manual(values= c(fillcolor))+
    labs(title= "Protein coverage", x="Counts", y= "Occurence")+
    geom_text(aes(label=counts), vjust=-0.2, size=4)+
    scale_x_continuous(breaks = 4:15 )
  
  return(bar)
}

#BOXPLOT
boxplot <- function (data,
                     width_boxplot=0.5,
                     color="black"){
  box <- ggplot2::ggplot(data, aes(x=label, y=intensity, fill=condition))+
    geom_boxplot(width=width_boxplot, color=color)+
    scale_fill_manual(values=color_palette)+
    theme_cuc()+
    theme(axis.text.x = element_text(angle = 30))+
    labs(title="Normalized data distribution", x="Label", y="Intensity")
  
  return(box)
}

#DENSITYPLOT
densityplot <- function (data,
                         color2="#481567FF",
                         color1="#21908CFF"){
  den <- data %>% 
    dplyr::group_by(gene_names) %>%
    dplyr::summarise(mean = mean(intensity, na.rm = TRUE),
                     missval = any(is.na(intensity))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(missing_value = dplyr::if_else(missval, "Missing", "Valid")) %>%
    dplyr::mutate(missing_value = factor(missing_value, levels = c("Valid", "Missing"))) %>%
    dplyr::group_by(missing_value) %>%
    
    ggplot2::ggplot(aes(x=mean, color=missing_value))+
    geom_density(linewidth=0.8, linetype='solid')+
    scale_color_manual(values = c(color1, color2))+
    theme_cuc()+
    labs(title="Density plot", x="Log2Intensity", y="Density")
  
  return(den)
}

#BARPLOT MISSING VALUES
barplot_missval <- function (data,
                             color="black",
                             color1="#481567FF",
                             color2="#21908CFF",
                             width_barplot=0.8){
  bar <- data %>% 
    dplyr::group_by(label) %>% 
    dplyr::mutate(bin_intensity = dplyr::if_else(bin_intensity == 1, "Valid", "Missing")) %>%
    dplyr::count(bin_intensity) %>% 
    dplyr::mutate(bin_intensity = as.factor(bin_intensity)) %>% 
    dplyr::rename(Data= bin_intensity) %>% 
    
    
    ggplot2::ggplot(aes(x=label, y=n, fill=Data))+
    geom_bar(stat="identity", width = width_barplot, color=color)+
    scale_fill_manual(values = c(color1, color2 ))+
    theme_cuc()+
    theme(axis.text.x = element_text(angle = 30))+
    labs(title= "Plot valid and missing data", x="Samples", y= "Counts")+
    geom_text(aes(label=n), position= position_stack(vjust = 0.5), size=3, show.legend = FALSE, color="white")
  
  return(bar)
}


#DENSITYPLOT IMPUTATION
den_imput <- function (data,
                       color1="#43BF71FF",
                       color2="#21908CFF",
                       color3="#35608DFF"){
  den <- data %>%
    dplyr::group_by(gene_names, condition) %>%
    dplyr::mutate(for_mean_imp = dplyr::if_else((sum(bin_intensity) / dplyr::n()) >= 0.75, TRUE, FALSE)) %>%
    dplyr::mutate(mean_grp = mean(intensity, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(imp_intensity = dplyr::case_when(
      bin_intensity == 0 & for_mean_imp ~ mean_grp,
      TRUE ~ as.numeric(intensity))) %>%
    dplyr::mutate(intensity = imp_intensity) %>% 
    dplyr::select(-c(for_mean_imp, mean_grp, imp_intensity))%>%
    dplyr::group_by(label) %>%
    # Define statistic to generate the random distribution relative to sample
    dplyr::mutate(mean = mean(intensity, na.rm = TRUE),
                  sd = sd(intensity, na.rm = TRUE),
                  n = sum(!is.na(intensity)),
                  total = nrow(data) - n) %>%
    dplyr::ungroup() %>%
    # Impute missing values by random draws from a distribution
    # which is left-shifted by parameter 'shift' * sd and scaled by parameter 'scale' * sd.
    dplyr::mutate(imp_intensity = dplyr::case_when(is.na(intensity) ~ rnorm(total,
                                                                            mean = (mean - 1.8 * sd), 
                                                                            sd = sd * 0.3),
                                                   TRUE ~ intensity)) %>%
    dplyr::mutate(intensity = imp_intensity) %>%
    dplyr::select(-c(mean, sd, n, total, imp_intensity)) %>%
    dplyr::group_by(condition) %>% 
    
    ggplot2::ggplot(aes(x=intensity, color=condition))+
    geom_density(linewidth=0.8, linetype='solid')+
    scale_color_manual(values = c(color1, color2, color3 ))+
    theme_cuc()+
    labs(title="Imputation plot", x="Log2Intensity", y="Density")
  
  return(den)
}

#HEATMAP
htmap <- function(data){
  map <- data %>%
    dplyr::select(gene_names, label, intensity) %>%
    tidyr::pivot_wider(names_from = label, values_from = intensity) %>%
    dplyr::filter(dplyr::if_all(.cols = dplyr::everything(), .fns = ~ !is.na(.x))) %>%
    tibble::column_to_rownames("gene_names") %>%
    cor() %>%
    round(digits = 2) %>% 
    Heatmap(name="Correlation",
            col = color_palette,
            column_title = "Correlation plot",
            column_title_gp=gpar(fontsize=18),
            cluster_rows = FALSE,
            cluster_columns = FALSE)
  
  return(map)
}

#PCA
#Make matrix
mat <- function(data) {
  mat <- data %>% 
    dplyr::select(gene_names, label, intensity) %>%
    tidyr::pivot_wider(id_cols = "gene_names",
                       names_from = "label",
                       values_from = "intensity") %>%
    tibble::column_to_rownames("gene_names") %>%
    as.matrix()
  return(mat)
}

#PCAPLOT
pca_plot <- function(data,
                     color1="#43BF71FF",
                     color2="#21908CFF",
                     color3="#35608DFF"){
  
  pca <- ggplot2::ggplot(data, aes(x=x, y=y),group=condition)+
  geom_point(size=3, shape=19, aes(color=condition))+
  theme_cuc()+
  geom_hline(yintercept = 0, linetype="longdash")+
  geom_vline(xintercept = 0, linetype="longdash")+
  scale_color_manual(values = c(color1, color2, color3))+
  labs(title="PCA", subtitle = "Principal component analysis", x="PC1", y="PC2")+
  geom_text(aes(label=replicate), size=3, position = "dodge", hjust=1.5)
  
  return(pca)
}


#DATASET_SIGNIFICANT
significant <- function(data, test){
  data <- data %>% 
    dplyr::select(gene_names, starts_with(test)) %>% 
    rename_at(vars(matches(test)), ~ str_remove(., paste0(test, "_"))) %>% 
    dplyr::filter(., !significant==FALSE) %>% 
    dplyr::mutate(regulation = dplyr::if_else(fold_change > 0, "Up", "Down")) %>% 
    dplyr::select(gene_names, fold_change, p_val, p_adj, regulation) %>%
    dplyr::mutate(gene_names = stringr::str_extract(gene_names, "[^;_]*")) %>% 
    dplyr::distinct(gene_names, .keep_all = TRUE) 
  return(data)
}


#FILTER_UP_DOWN_REGULATED
sig_up_down <- function(data, remove){
  data <- data %>% 
    dplyr::filter(., !regulation == remove)
  return(data)
}


#NODES STRING (diventa inutile perchè quando andiamo a cercare gli edges la ricerca viene già fatta tramite gene_name.
#Di conseguenza basta modificare la sig table per ottenere la lista di nodi)
nodes <- function(data, nodesize=5){
  data <- data %>% 
    dplyr::pull(gene_names) %>% 
    rba_string_map_ids(species = 9606, echo_query = TRUE) %>%
    dplyr::select(stringId, preferredName, annotation ) %>% 
    dplyr::rename(name=preferredName) %>% 
    dplyr::left_join(sig, by= c("name"="gene_names")) %>% 
    dplyr::select(-c(stringId, annotation, fold_change)) %>% 
    dplyr::rename(value =p_adj, size = p_val, grp = regulation) %>% 
    dplyr::mutate(value = -log10(value)) %>% 
    dplyr::mutate(size = -log10(size)*nodesize) %>% 
    dplyr::mutate(grp = as.factor(grp)) %>% 
    tidyr::as_tibble()
  
  return(data)
}

#EDGES STRING
edges <- function(data, score=0.4, edgesize=5){
  data <- data %>% 
    dplyr::pull(name) %>%
    rba_string_interactions_network(species = 9606) %>% 
    dplyr::mutate(Names= paste0(preferredName_A, "_", preferredName_B)) %>% 
    dplyr::distinct(Names, .keep_all = TRUE) %>% 
    dplyr::rename(FROM = preferredName_A, TO=preferredName_B) %>%
    dplyr::mutate(totscore = escore+dscore) %>% 
    dplyr::filter(!totscore== "0") %>%
    dplyr::select(FROM,TO, escore, dscore) %>%
    #se facessi tidyr::unite() potrei unire due colonne in una unica e poi fare il distinct
    dplyr::mutate(Score1 = (escore-0.041)*(1-0.041)) %>% 
    dplyr::mutate(Score2 = (dscore-0.041)*(1-0.041)) %>% 
    dplyr::mutate(Score_combin = 1-(1-Score1)*(1-Score2)) %>% 
    dplyr::mutate(tot_score= Score_combin+0.041*(1-Score_combin)) %>% 
    dplyr::filter(tot_score>score) %>% #Filtro sullo score
    dplyr::mutate(tot_score = tot_score*edgesize) %>% 
    dplyr::select(FROM,TO, tot_score) %>% 
    dplyr::rename(source=FROM, target=TO, value=tot_score)
  
  return(data)
}


#COMPLEXES CORUM
complexes_corum <- function(data = up_down) {
  all_complexes <-
    import_omnipath_complexes(resources = "CORUM") #Tutte le interactions di CORUM
  query <- data$gene_names #Dataset di partenza
  my_complexes <- unique(get_complex_genes(all_complexes, query,
                                           total_match = FALSE))
  return(my_complexes)
}

#NODES CORUM
nodes_corum <- function(data, nodesize){
  data <- data %>% 
    dplyr::select(components_genesymbols) %>% 
    tidyr::separate_rows(components_genesymbols, sep = "_") %>%
    dplyr::filter(components_genesymbols %in% query) %>% 
    dplyr::distinct(components_genesymbols, .keep_all = TRUE) %>% #Tolgo nomi duplicati
    dplyr::rename(name=components_genesymbols) %>% 
    dplyr::left_join(up_down, by=c("name"="gene_names")) %>%
    dplyr::rename(value = p_adj, size = p_val, grp = regulation) %>% 
    dplyr::mutate(grp = "corum") %>%
    dplyr::mutate(value = -log10(value)) %>%
    dplyr::mutate(size = -log10(size) * nodesize) %>%  
    dplyr::select(name, value, size, grp) %>% 
    dplyr::distinct(name, .keep_all = TRUE) 
  
  return(data)
}

#Quando voglio far capire dove mettere il termine della pipe precedente, mi basta mettere il .
#Altrimenti me lo mette come primo termine

#FILTER_NODES_DEGREE_ZERO
filter_nodes <- function(datanodes, data_edges){
  data <- datanodes %>% 
    left_join(data_edges, by=c("name"="source")) %>% 
    left_join(data_edges, by=c("name"="target")) %>% 
    dplyr::mutate(source = dplyr::if_else(is.na(source), 0, 1)) %>% 
    dplyr::mutate(target = dplyr::if_else(is.na(target), 0, 1)) %>% 
    dplyr::mutate(Exists = source+target) %>% 
    dplyr::filter(!Exists== "0") %>% 
    dplyr::select(name, size, value, grp) %>% 
    dplyr::distinct(name, .keep_all = TRUE)
  
  return(data)
}

##TUTTI I NODI DA RAPPRESENTARE
all_nodes <- function(nodes_string, nodes_corum){
  data <- nodes_string %>% 
    dplyr::mutate(grp="string") %>% 
    dplyr::bind_rows(nodes_corum) %>% as.data.frame() 
  
  return(data)
}

##EDGES CORUM
edges_corum <- function(data_filtered, data, data_edges) {
  query <- data_filtered$gene_names #Dataset di partenza
  
  data <- data %>%
    dplyr::select(name, components_genesymbols) %>%
    tidyr::separate_rows(components_genesymbols, sep = "_") %>%
    dplyr::filter(components_genesymbols %in% query) %>%
    dplyr::rename(source = components_genesymbols) %>%
    unique() %>%  #In questo modo seleziono solo una riga
    get_dupes(name) %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(source) %>%
    dplyr::filter(count == max(count)) %>%
    dplyr::distinct(source, .keep_all = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(complex = name) %>%
    dplyr::group_by(name) %>%
    dplyr::right_join(data_edges) %>%
    dplyr::select(source, target, value, name) %>% 
    dplyr::mutate(value=10)
  
  return(data)
  
}



#PLOT_NETWORK
plot_net <- function(datanode, dataedge, animation=FALSE, layout="force"){
  p <- echarts4r::e_charts(animation=FALSE) %>% 
    echarts4r::e_graph(roam=TRUE, 
                       layout = layout,
                       force= list( initLayout = "circular", 
                                    repulsion=800,
                                    edgeLength=50,
                                    layoutAnimation=animation),
                       itemStyle=list(opacity=0.9),
                       autoCurveness = TRUE,
                       emphasis=list(focus="adjacency")) %>% 
    echarts4r::e_graph_nodes(nodes=datanode, names = name,  value=value, size=size, category=grp) %>% 
    echarts4r::e_graph_edges(edges=dataedge, source=source, target=target, value=value, size = value) %>%
    echarts4r::e_labels(datanode$name, font_size=4)%>%
    echarts4r:: e_title("Network", "Up & Down Regulated") %>%
    echarts4r::e_tooltip()
  
  return(p)
}

#PLOT_ALL
plot_all <- function(data, test, remove, score=0.4, animation=FALSE, layout="force", no_edge=TRUE){
  sig <- significant(data, test)
  up_down <- sig_up_down(data=sig, remove)
  nodes <- nodes(data=up_down)
  edges <- edges(data= nodes, score)
  if (no_edge) {nodes<- filter_nodes (datanodes =nodes, data_edges = edges)}
  p <- plot_net(datanode=nodes, dataedge=edges, animation, layout)
  
  return(p)
}

##CORUM EDGES COLOR
color_edge <- function(list, edges) {
  n_edges <- nrow(edges)
  for (i in 1:n_edges) {
    source <-
      list %>%
      purrr::pluck("x", "opts", "series", 1, "links", i, "source")
    
    target <-
      list %>% purrr::pluck("x", "opts", "series", 1, "links", i, "target")
    
    width <-
      list %>% purrr::pluck("x", "opts", "series", 1, "links", i, "lineStyle", "width")
    
    val <- as.numeric(width) %>% round(0)
    
    color <-
      edges %>%
      dplyr::filter(source == source, target == target, value == val) %>%
      pull(color)
    
    list <-
      purrr::modify_in(list,
                       list("x", "opts", "series", 1, "links", i, "lineStyle"),
                       \(x) c(x, c(color = color)))
  }
  return(list)
  
}


### OTHER DATABASES SEARCH

## Intact
intact_edge <- function(gene_vector) {
  intact_list <- rba_reactome_interactors_psicquic(genes, resource = "IntAct", details = TRUE)
  
  all_interactions <- purrr::map(.x = 1:length(gene_vector),
                                 .f = ~ intact_list$entities[[.x]])
  
  res <- purrr::map(.x = all_interactions, .f = ~ keep_only_if_present(.x)) %>%
    purrr::compact() %>%
    purrr::reduce(bind_rows) %>%
    dplyr::select(source, alias, score, accURL, evidences) %>%
    tidyr::drop_na() %>%
    dplyr::filter(source != alias)
  
  return(res)
}


#Keep only nodes present in my dataset
keep_only_if_present <- function(lista) {
  n_interaction <- lista$count
  name_source <- lista$acc
  
  res <-
    purrr::map(
      .x = 1:n_interaction,
      .f = ~ lista$interactors[[.x]] %>% purrr::list_modify(source = name_source)
    ) %>%
    purrr::keep(~ all(.x$alias %in% pull(nodes, name)))
  
  return(res)
}

