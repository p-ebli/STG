# 
#   STG diagram using getTransitionTable()
# 
# 
# 



library(BoolNet)
library(igraph)

# Set working directory 
setwd("Path\\to\\your\\directory")


extract_boolean_expressions <- function(file_name) {
  lines <- readLines(file_name, warn = FALSE)
  
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  return(lines)
}

# List files and sort them in natural numeric order.
txt_files <- list.files(pattern = "\\.txt$")
# This sorts based on the numeric part (i.e. numbers in file names)
txt_files <- txt_files[order(as.numeric(gsub("\\D", "", txt_files)))]

all_results <- list()

for (file_name in txt_files) {
  cat("-------------------------------------------------\n")
  cat("Processing network:", file_name, "\n")
  
  
  net <- loadNetwork(file_name)
  print(net)
  
  
  boolean_expressions <- extract_boolean_expressions(file_name)
  cat("Extracted Boolean expressions:\n")
  print(boolean_expressions)
  

  if (length(boolean_expressions) > 1) {
    assignments <- boolean_expressions[-1]  # remove header
    all_zero <- TRUE
    for (line in assignments) {
      tokens <- strsplit(line, ",")[[1]]
      if (length(tokens) >= 2) {
        value <- trimws(tokens[2])
        if (value != "0") {
          all_zero <- FALSE
          break
        }
      }
    }
    if (all_zero) {
      cat("Skipping network", file_name, "- all genes are set to 0.\n\n")
      next
    }
  }
  

  att <- getAttractors(net, type = "synchronous", returnTable = TRUE)
  
 
  tt_text <- capture.output(getTransitionTable(att))
  
  
  transition_lines <- grep("=>", tt_text, value = TRUE)
  
 
  edges_list <- lapply(transition_lines, function(line) {
    parts <- strsplit(line, "=>")[[1]]
    if (length(parts) < 2) return(NULL)
    from_state <- trimws(parts[1])
    tokens <- strsplit(trimws(parts[2]), "\\s+")[[1]]
    if (length(tokens) < 1) return(NULL)
    to_state <- tokens[1]
    return(c(from = from_state, to = to_state))
  })
  

  edges_list <- Filter(Negate(is.null), edges_list)
  if (length(edges_list) == 0) {
    cat("No transitions found for network:", file_name, "\n")
    next
  }
  edges_mat <- do.call(rbind, edges_list)
  edges_df <- as.data.frame(edges_mat, stringsAsFactors = FALSE)
  edges_df <- edges_df[edges_df$from != "" & edges_df$to != "", ]
  
  cat("\nConstructed edges data frame:\n")
  print(edges_df)
  
  g <- graph_from_data_frame(edges_df, directed = TRUE)
  
  
  all_results[[file_name]] <- list(
    network = net,
    attractors = att,
    boolean_expressions = boolean_expressions,
    edges_df = edges_df,
    graph = g
  )
  
  cat("Finished processing network:", file_name, "\n\n")
}


ordered_names <- names(all_results)[order(as.numeric(gsub("\\D", "", names(all_results))))]
all_results <- all_results[ordered_names]

# ------------------------------------------
# PDF generation
# ------------------------------------------
pdf("State_Transition_Graphs.pdf", onefile = TRUE, width = 12, height = 8)
for (nm in names(all_results)) {
  
  par(mfrow = c(1, 2), mar = c(4, 4, 4, 2))
  
  
  g <- all_results[[nm]]$graph
  plot(g,
       vertex.shape = "circle",
       vertex.size = 40,
       vertex.color = "lightblue",
       vertex.label.color = "black",
       edge.arrow.size = 0.5,
       main = paste("State Transition Mapping:", nm))
  
  
  plot.new()
  boolean_expressions <- all_results[[nm]]$boolean_expressions
  if (length(boolean_expressions) > 0) {
    n_lines <- length(boolean_expressions)
   
    y_coords <- 0.95 - (0:(n_lines - 1)) * 0.05
    for (i in seq_along(boolean_expressions)) {
      text(0, y_coords[i], labels = boolean_expressions[i], adj = c(0, 1), cex = 0.8)
    }
  } else {
    text(0, 0.9, "No Boolean expressions found.", adj = c(0, 1), cex = 0.8)
  }
}
dev.off()

cat("All graphs along with Boolean expressions have been saved to 'State_Transition_Graphs.pdf'\n")
