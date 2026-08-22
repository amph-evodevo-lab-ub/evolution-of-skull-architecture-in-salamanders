library (ape)
library (phytools)
library(geiger)
library(dplyr)
library(igraph)
library(scales)


#######################################Data preparation#############################
###Read NEXUS tree file
tree <- read.nexus("Tree_salamandridae_cranial_articulations.nex")
print (tree)

###Read character states and mean values for each species
character_states = read.csv("R_data_salamandridae_cranial_articulations.csv", header = T)
character_states


###Prune phylogenetic tree only to contain species with character states data
traits <- names(character_states)[
  sapply(character_states, is.numeric)
]
traits

species_to_keep <- character_states$species
species_to_keep

pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, species_to_keep))
pruned_tree
plot(pruned_tree)


#####################Mapping character states onto phylogeny with ancestral states###########
###Prepare results table

traits <- setdiff(colnames(character_states), "species")


anc_results <- data.frame(
  trait = character(),
  root_ace_BM = numeric(),
  root_ace_OU = numeric(),
  stringsAsFactors = FALSE
)


###Make a loop to analyze all traits 
n_traits <- length(traits)

for (tr in traits) {
  
  message("Processing trait: ", tr)
  
 
  trait_vector <- setNames(
    character_states[[tr]],
    character_states$species
  )
  
  trait_vector <- trait_vector[pruned_tree$tip.label]
  trait_vector <- trait_vector[!is.na(trait_vector)]
  
  tree_tr <- drop.tip(
    pruned_tree,
    setdiff(pruned_tree$tip.label, names(trait_vector))
  )
  
  root_node <- Ntip(tree_tr) + 1
  

  anc_bm <- fastAnc(
    tree_tr,
    trait_vector
  )
  
  root_bm <- anc_bm[as.character(root_node)]
  
  
  anc_ou <- anc.ML(
    tree_tr,
    trait_vector,
    model = "OU"
  )
  
  root_ou <- anc_ou$ace[as.character(root_node)]
  
 
  anc_results <- rbind(
    anc_results,
    data.frame(
      trait = tr,
      root_ace_BM = root_bm,
      root_ace_OU = root_ou
    )
  )
  

  cm <- contMap(
    tree_tr,
    trait_vector,
    plot = FALSE,
    res = 600
  )

  tiff(paste0("contMap_", tr, ".tiff"),
       width = 6,
       height = 8,
       units = "in",
       res = 600)
  
  plot(
    cm,
    fsize = 0.8,
    lwd = 2,
    
    main = tr
  )
  
  nodelabels(
    node = root_node,
    text = round(root_ou, 2),
    frame = "circle",
    bg = "white",
    cex = 0.6
  )
  
  dev.off()
}
  

write.csv(
  anc_results,
  "Ancestral_states_root_BM_OU.csv",
  row.names = FALSE
)



######################Calculate phylogenetic signal################

###Prepare results table

results_signal <- data.frame(
  trait = traits,
  lambda = NA,
  lambda_p = NA,
  K = NA,
  K_p = NA
)



###Make a loop to analyze all traits 

for (i in seq_along(traits)) {
  trait <- traits[i]
  
  vec <- setNames(character_states[[trait]], character_states$species)
  

  vec <- vec[!is.na(vec)]

  vec <- vec[pruned_tree$tip.label]

  res_lambda <- phylosig(pruned_tree, vec, method = "lambda", test = TRUE)
  
 
  res_K <- phylosig(pruned_tree, vec, method = "K", test = TRUE)
  
 
  results_signal$lambda[i]   <- res_lambda$lambda
  results_signal$lambda_p[i] <- res_lambda$P
  
  results_signal$K[i]   <- res_K$K
  results_signal$K_p[i] <- res_K$P
}

# FDR correction across all p-values
results_signal <- results_signal %>%
  mutate(
    lambda_p_FDR = p.adjust(lambda_p, method = "fdr"),
    K_p_FDR = p.adjust(K_p, method = "fdr")
  )


results_signal
results_signal <- results_signal |> 
  mutate(across(where(is.numeric), ~ round(.x, 4)))
write.csv(results_signal, "phylogenetic_signal_results_corr", row.names = FALSE)


###############################Evolutionary model check######################  
results <- list()

n_traits <- length(traits)
n_traits

###Make a loop to analyze all traits 

for (i in seq_along(traits)) {
  
  tr <- traits[i]
  message("[", i, "/", n_traits, "] Storing results for trait: ", tr)
  
  
  trait_vector <- setNames(
    character_states[[tr]],
    character_states$species
  )
  
  trait_vector <- trait_vector[pruned_tree$tip.label]
  
  stopifnot(
    length(trait_vector) == length(pruned_tree$tip.label),
    !any(is.na(trait_vector))
  )

  fitBM <- fitContinuous(pruned_tree, trait_vector, model = "BM")
  fitOU <- fitContinuous(pruned_tree, trait_vector, model = "OU")
  fitEB <- fitContinuous(pruned_tree, trait_vector, model = "EB")
  
  AICc <- c(
    BM = fitBM$opt$aicc,
    OU = fitOU$opt$aicc,
    EB = fitEB$opt$aicc
  )
  
  delta   <- AICc - min(AICc)
  weights <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))
  
  
  anc_ou <- anc.ML(
    pruned_tree,
    trait_vector,
    model = "OU"
  )
  
  root_node  <- Ntip(pruned_tree) + 1
  root_state <- anc_ou$ace[as.character(root_node)]
  

  results[[tr]] <- list(
    root_state = root_state,
    theta      = fitOU$opt$theta,
    AICc       = AICc,
    weights    = weights
  )
}

summary_table <- do.call(rbind, lapply(names(results), function(tr) {
  r <- results[[tr]]
  
 
  
  data.frame(
    trait = tr,
    root_state = r$root_state,
    AICc_BM = r$AICc["BM"],
    AICc_OU = r$AICc["OU"],
    AICc_EB = r$AICc["EB"],
    weight_BM = r$weights["BM"],
    weight_OU = r$weights["OU"],
    weight_EB = r$weights["EB"],
    stringsAsFactors = FALSE
  )
}))
row.names(summary_table) <- NULL

summary_table
write.csv(summary_table, "model_check.csv", row.names = FALSE)


############################PIC analysis#########################

###data preparation
df <- character_states %>%
  filter(species %in% pruned_tree$tip.label) %>%
  distinct(species, .keep_all = TRUE) %>%
  arrange(match(species, pruned_tree$tip.label))

stopifnot(all(df$species == pruned_tree$tip.label))

trait_pairs <- combn(traits, 2, simplify = FALSE)

###prepare results table 
pic_results <- data.frame(
  trait_x = character(),
  trait_y = character(),
  slope   = numeric(),
  p_value = numeric(),
  r2      = numeric(),
  stringsAsFactors = FALSE
)

###make a loop to analyze all trait pairs

for (pair in trait_pairs) {
  
  tr1 <- pair[1]
  tr2 <- pair[2]
  
  message("PIC: ", tr1, " vs ", tr2)
  
  x <- df[[tr1]]
  y <- df[[tr2]]
  
  
  if (any(is.na(x)) || any(is.na(y))) {
    message("  Skipped (NA values)")
    next
  }
  
  
  pic_x <- try(pic(x, pruned_tree), silent = TRUE)
  pic_y <- try(pic(y, pruned_tree), silent = TRUE)
  
  if (inherits(pic_x, "try-error") || inherits(pic_y, "try-error")) {
    message("  PIC failed")
    next
  }
  
  
  model <- lm(pic_y ~ pic_x - 1)
  sm <- summary(model)
  
  
  pic_results <- rbind(
    pic_results,
    data.frame(
      trait_x = tr1,
      trait_y = tr2,
      slope   = sm$coefficients[1, "Estimate"],
      p_value = sm$coefficients[1, "Pr(>|t|)"],
      r2      = sm$r.squared,
      stringsAsFactors = FALSE
    )
  )
}
pic_results
pic_results$p_adj_BH <- p.adjust(pic_results$p_value, method = "BH")
pic_results <- pic_results |> 
  mutate(across(where(is.numeric), ~ round(.x, 4)))
write.csv(pic_results, "PIC_trait_pairs-1jun26.csv", row.names = FALSE)




###############Network of statistically significant correlations among examined cranial articulation traits
#####################based on phylogenetically independent contrasts#########################

###exclude traits other than 10 analysed characters
trait_names <- setdiff(
  colnames(character_states),
  c("species", "CS(mm)", "corrected_thickness", "CAI", "aquatic_period")
)

###PIC for 10 characters
pairwise_results <- data.frame(
  from = character(),
  to   = character(),
  p    = numeric(),
  R2   = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:(length(trait_names) - 1)) {
  for (j in (i + 1):length(trait_names)) {
    
    tr1 <- trait_names[i]
    tr2 <- trait_names[j]
    
    v1 <- character_states[[tr1]]
    v2 <- character_states[[tr2]]
    
    names(v1) <- names(v2) <- character_states$species
    
    v1 <- v1[pruned_tree$tip.label]
    v2 <- v2[pruned_tree$tip.label]
    
    
    keep <- complete.cases(v1, v2)
    v1 <- v1[keep]
    v2 <- v2[keep]
    tree_sub <- drop.tip(pruned_tree, pruned_tree$tip.label[!keep])
    
    
    if (length(v1) < 5 || sd(v1) == 0 || sd(v2) == 0) next
    
    pic1 <- pic(v1, tree_sub)
    pic2 <- pic(v2, tree_sub)
    
    fit <- lm(pic1 ~ pic2 - 1)
    
    pairwise_results <- rbind(
      pairwise_results,
      data.frame(
        from = tr1,
        to   = tr2,
        p    = summary(fit)$coefficients[1,4],
        R2   = summary(fit)$r.squared
      )
    )
  }
}

pairwise_results$p_adj <- p.adjust(pairwise_results$p, method = "BH")

sig_pairs <- subset(pairwise_results, p_adj < 0.05)

####diagram graph
nodes <- data.frame(
  name = unique(c(sig_pairs$from, sig_pairs$to)),
  stringsAsFactors = FALSE
)

edges <- rbind(
  data.frame(
    from = sig_pairs$from,
    to   = sig_pairs$to,
    weight = sig_pairs$R2,
    p = sig_pairs$p_adj
  )
)

g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

node_colors <- c(
  "premaxilla_maxilla" = "red",
  "nasal_frontal" = "red",
  "prefrontal_frontal" = "red",
  
  "frontal_parietal" = "cadetblue2",
  "skull_roof_left_right" = "cadetblue2",
  
  "parietal_otic-occipital" = "green",
  "orbitosphenoid_connections" = "green",
  
  "vomer_parasphenoid" = "royalblue3",
  "parasphenoid_otic-occipital" = "royalblue3",

  "squamosum_otic-occipital" = "royalblue3"
)

trait_labels <- c(
  "premaxilla_maxilla" = "1",
  "nasal_frontal" = "2",
  "prefrontal_frontal" = "3",
  "frontal_parietal" = "4",
  "skull_roof_left_right" = "5",
  "parietal_otic-occipital" = "6",
  "orbitosphenoid_connections" = "7",
  "vomer_parasphenoid" = "8",
  "parasphenoid_otic-occipital" = "9",
  "squamosum_otic-occipital" = "10"
)


if (!"squamosum_otic-occipital" %in% V(g)$name) {
  g <- add_vertices(
    g,
    1,
    name = "squamosum_otic-occipital"
  )
}

###node settings###
V(g)$color <- node_colors[V(g)$name]
V(g)$size <- 18
V(g)$frame.color <- "black"
V(g)$label <- trait_labels[V(g)$name]
V(g)$label.color <- ifelse(
  V(g)$color == "royalblue3",
  "white",
  "black"
)

set.seed(123)
lay <- layout_with_fr(g)


id10 <- which(V(g)$name == "squamosum_otic-occipital")
id8  <- which(V(g)$name == "vomer_parasphenoid")
id9  <- which(V(g)$name == "parasphenoid_otic-occipital")

center_x <- mean(c(lay[id8,1], lay[id9,1]))
center_y <- mean(c(lay[id8,2], lay[id9,2]))


lay[id10,1] <- center_x - 0.8
lay[id10,2] <- center_y - 0.3


###edge settings###
E(g)$width <- scales::rescale(E(g)$weight, to = c(1, 8))
E(g)$color <- "gray45"



###plot###

par(mar = c(2,2,2,2))

plot(
  g,
  layout = lay,
  vertex.label.cex = 1
)

#legend#
r2_vals <- E(g)$weight

r2_min <- min(r2_vals)
r2_max <- max(r2_vals)

scale_width <- function(x) {
  scales::rescale(x, to = c(1, 8), from = range(r2_vals))
}

lwd_vals <- c(
  scale_width(r2_min),
  scale_width(r2_max)
)

legend(
  "topleft",
  legend = c(
    "Strength",
    paste0("R² min = ", round(r2_min, 2)),
    paste0("R² max = ", round(r2_max, 2))
  ),
  lwd = c(NA, lwd_vals[1], lwd_vals[2]),
  col = c(NA, "gray45", "gray45"),
  bty = "n",
  cex = 0.7
)

