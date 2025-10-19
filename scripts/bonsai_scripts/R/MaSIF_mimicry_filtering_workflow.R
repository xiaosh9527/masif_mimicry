library(dplyr)

# R script to filter MaSIF mimicry truncated binders by various postprocessing
# metrics, and save to a filtered .csv of truncated binders. 

################################################################################

# ----- 1. Read input dataframes ------

# df containing the unfiltered truncation scores
df_trunc_score <- read.csv("./results/proc_trunc_86_6H0F_C_B.csv")

# df containing the raw MaSIF mimicry scores
df_mimicry <- read.csv("./results/6h0f_C_B_hp_afdb_merizo_split_search_v2_masif_score.csv")

# df containing the filtering thresholds
df_filter <- read.csv("./results/trunc_filters_mimicry_6H0F_C_B_250325.csv")


# ----- 2. Calculate additional postprocessing metrics ---------

### iface-to-terminus metrics:
# Function to calculate iface-to-terminus metrics for a single row
iface_to_term_metrics <- function(iface_resi, resi_start, resi_end) {
  # Initialize return values
  min_dist <- NA
  mean_dist <- NA
  
  # Ensure that the iface residue string is not missing or empty
  if (!is.na(iface_resi) && iface_resi != "") {
    # Split the underscore-delimited string and convert to numeric
    iface_nums <- as.numeric(unlist(strsplit(as.character(iface_resi), "_")))
    
    # Calculate the distances for each iface residue to the N- and C-terminus
    dist_to_n <- iface_nums - resi_start       # Distance from N-terminus
    dist_to_c <- resi_end - iface_nums         # Distance from C-terminus
    
    # For each residue, take the smaller distance
    distances <- pmin(dist_to_n, dist_to_c)
    
    # The minimum iface distance to a terminus is the smallest of these values
    min_dist <- min(distances)
    
    # Calculate the mean residue number of the iface residues
    mean_iface <- mean(iface_nums)
    
    # Compute the distance of the mean to each terminus
    mean_dist_to_n <- mean_iface - resi_start
    mean_dist_to_c <- resi_end - mean_iface
    
    # Record the smaller of these two distances
    mean_dist <- min(mean_dist_to_n, mean_dist_to_c)
  }
  
  return(list(min_dist = min_dist, mean_dist = mean_dist))
}

# Apply the function to df_trunc_score using mapply
metrics_results <- mapply(iface_to_term_metrics, 
                         df_trunc_score$trunc_iface_residues,
                         df_trunc_score$trunc_resi_start, 
                         df_trunc_score$trunc_resi_end, 
                         SIMPLIFY = FALSE)

# Extract the results and add to dataframe
df_trunc_score$min_iface_dist_terminus <- sapply(metrics_results, function(x) x$min_dist)
df_trunc_score$mean_iface_dist_terminus <- sapply(metrics_results, function(x) x$mean_dist)

### uniprot_exclude:
# (Exclude uniprot ids that frequently appear and makes no biological sense. uniprot to be excluded will have a value of 1)
df_uniprot_to_exclude <- read.csv("./uniprot_to_exclude.csv")

### domain_exclude:
# (Exclude domains that frequently appear and makes no biological sense, but not currently in use, so df_domain_exclude is empty)
df_domains_to_exclude <- read.csv("./domains_to_exclude.csv")

# remove old uniprot_exclude columns if they exist
df_trunc_score <- df_trunc_score %>%
  select(-any_of(c( "uniprot_exclude", "domain_exclude")))

# Merge these tables
df_trunc_score <- merge(
  df_trunc_score, 
  df_uniprot_to_exclude,
  by = "uniprot_id", 
  all.x = TRUE
)
df_trunc_score <- merge(
  df_trunc_score, 
  df_domains_to_exclude,
  by = "matched_protein", 
  all.x = TRUE
)
df_trunc_score <- 
  mutate(df_trunc_score,
         uniprot_exclude = if_else(is.na(uniprot_exclude), 0, uniprot_exclude),
         domain_exclude = if_else(is.na(domain_exclude), 0, domain_exclude)
  )



# ----- 3. Merge df_trunc_score with df_mimicry to get MaSIF metrics ------

# Add MaSIF scores to df_trunc_score by matched_protein 
df_mimicry$matched_protein <- gsub("_A", "", df_mimicry$P1_id)
df_mimicry <- select(df_mimicry,
                    matched_protein,
                    P1_source_TMscore, P1_source_iface, 
                    P2_source_TMscore, P2_source_iface, 
                    MaSIF.score, aa_match, ss_match)
df_mimicry_unique <- df_mimicry %>%
  group_by(matched_protein) %>%
  slice_max(MaSIF.score, n = 1, with_ties = FALSE) %>%
  ungroup()

# merge dataframe
df_trunc_score_merged <- merge(df_trunc_score, df_mimicry_unique,
                               by = "matched_protein",
                               all.x = TRUE)


# ----- 4. Filter df_trunc_score_merged by the filtering thresholds in df_filters ------

# Create a function to apply filtering thresholds
apply_filters <- function(df, filter_df) {
  filtered_df <- df
  
  # Loop through each filter threshold
  for (i in 1:nrow(filter_df)) {
    metric <- filter_df$metric[i]
    min_val <- filter_df$min[i]
    max_val <- filter_df$max[i]
    
    # Check if the metric column exists in the dataframe
    if (metric %in% colnames(filtered_df)) {
      # Apply the filter (keep rows where metric is between min and max)
      filtered_df <- filtered_df[filtered_df[[metric]] >= min_val & 
                                filtered_df[[metric]] <= max_val & 
                                !is.na(filtered_df[[metric]]), ]
      
      cat(sprintf("Filtered by %s: %d rows remaining\n", metric, nrow(filtered_df)))
    } else {
      cat(sprintf("Warning: Column %s not found in dataframe\n", metric))
    }
  }
  
  return(filtered_df)
}

# Apply all filters to the merged dataframe
df_filtered <- apply_filters(df_trunc_score_merged, df_filter)


# ----- 5. Save the filtered results ------

# Save the filtered dataframe
write.csv(df_filtered, 
          file = "./results/proc_trunc_86_6H0F_C_B_250325_filtered.csv",
          row.names = FALSE)














