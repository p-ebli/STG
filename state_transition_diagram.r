#############################################################################################
###    this code loads networks from a directory                                          ###
###    based on the number of genes in the network and their values generates states      ###
###      if in the network, a gene is set to 0, in all states has value 0                 ###
###      if all genes are set to 0, network is skipped                                    ###
###    computes attractors with getAttractors()                                           ###
###    computes path to attractor with getPathToAttractor using states and attractors     ###
###    generates state transition diagrams                                                ###
###    save as pdf file                                                                   ###
#############################################################################################



library(BoolNet)
library(igraph)

# Set working directory
setwd("Path\\to\\your\\directory")

# Convert a numeric state vector into a compact string.
vector_to_state <- function(state_vector) {
  paste(state_vector, collapse = "")
}

# Function to generate valid states based on fixed gene constraints.
generate_valid_states <- function(net, fixed_constraints) {
  genes <- net$genes
  possible_values <- lapply(genes, function(gene) {
    if (gene %in% names(fixed_constraints)) {
      fixed_constraints[[gene]]
    } else {
      c(0, 1)
    }
  })
  names(possible_values) <- genes
  states_df <- expand.grid(possible_values)
  states_df <- states_df[, genes, drop = FALSE]
  states <- lapply(seq_len(nrow(states_df)), function(i) as.numeric(states_df[i, ]))
  return(states)
}

# ------------------------------
# Process all .txt files in directory
# ------------------------------
all_results <- list()
txt_files <- list.files(pattern = "\\.txt$")

for (file_name in txt_files) {
  cat("Processing network:", file_name, "\n")
  
  # Load the network
  net <- loadNetwork(file_name)
  print(net)
  
  # Read network description for fixed constraints
  network_df <- read.csv(file_name, header = TRUE, stringsAsFactors = FALSE)
  fixed_constraints <- list()
  for (j in seq_len(nrow(network_df))) {
    gene <- trimws(network_df$targets[j])
    factor_val <- trimws(network_df$factors[j])
    if (factor_val %in% c("0", "1")) {
      fixed_constraints[[gene]] <- as.numeric(factor_val)
    }
  }
  
  cat("Fixed constraints for this network:\n")
  print(fixed_constraints)
  
  # Skip network if all genes fixed to 0
  if (length(fixed_constraints) == length(net$genes) &&
      all(unlist(fixed_constraints) == 0)) {
    cat("Skipping network", file_name, "All gene values = 0.\n\n")
    next
  }
  
  # Compute attractors
  att <- getAttractors(net, "synchronous", returnTable = TRUE)
  print(att)
  
  # Generate valid initial states
  valid_states <- generate_valid_states(net, fixed_constraints)
  cat("Valid states:\n")
  state_labels <- sapply(valid_states, vector_to_state)
  print(state_labels)
  
  # Compute path to attractors for each state
  paths <- list()
  state_mapping <- list()
  for (state in valid_states) {
    if (length(state) != length(net$genes)) stop("State length mismatch!")
    state_label <- vector_to_state(state)
    cat("Computing path for state:", state_label, "\n")
    path <- getPathToAttractor(net, state, "all")
    paths[[state_label]] <- path
    attractor_state <- vector_to_state(as.numeric(path[nrow(path), ]))
    state_mapping[[state_label]] <- attractor_state
    cat("Path for state", state_label, ":\n")
    print(path)
  }
  
  # Build edge list for igraph
  edge_list <- c()
  cat("State mapping:\n")
  for (init_state in names(state_mapping)) {
    attractor <- state_mapping[[init_state]]
    cat(init_state, "results in", attractor, "\n")
    edge_list <- c(edge_list, init_state, attractor)
  }
  
  g <- graph(edges = edge_list, directed = TRUE)
  
  # Optional interactive plot
  plot(g,
       vertex.shape = "circle",
       vertex.size = 40,
       vertex.color = "lightblue",
       vertex.label.color = "black",
       edge.arrow.size = 0.5,
       main = paste("State-to-Attractor Mapping for", file_name))
  
  all_results[[file_name]] <- list(
    network = net,
    attractors = att,
    paths = paths,
    state_mapping = state_mapping,
    graph = g,
    description = network_df
  )
  
  cat("Finished processing network:", file_name, "\n\n")
}

# ------------------------------
# Save all graphs to a PDF
# ------------------------------
pdf("State_to_Attractor_Graphs.pdf", onefile = TRUE, width = 12, height = 8)
for(nm in names(all_results)) {
  par(mfrow = c(1, 2), mar = c(4, 4, 4, 2))
  
  g <- all_results[[nm]]$graph
  plot(g,
       vertex.shape = "circle",
       vertex.size = 40,
       vertex.color = "lightblue",
       vertex.label.color = "black",
       edge.arrow.size = 0.5,
       main = paste("State-to-Attractor Mapping:\n", nm))
  
  plot.new()
  network_info <- all_results[[nm]]$description
  info_text <- paste(capture.output(print(network_info)), collapse = "\n")
  text(0, 1, info_text, adj = c(0, 1), cex = 0.8)
}
dev.off()
cat("All graphs have been saved to 'State_to_Attractor_Graphs.pdf'\n")
