library(GenomicRanges)

# Create a GRanges object representing your exons (replace with your actual data)
exon_ranges <- GRanges(
  seqnames = Rle("chr1", 5),    # Chromosome names for exons
  ranges = IRanges(
    start = c(100, 200, 300, 400, 500),  # Start positions of exons
    end = c(150, 250, 350, 450, 550)    # End positions of exons
  )
)

# Define the position of the random site (replace with your actual position)
site <- 250

# Calculate the distances from the site to ALL exons manually
start_distances_y <- abs(start(exon_ranges) - site)
stop_distances_y <- abs(end(exon_ranges) - site)

y_values <- c() 
for (i in 1:length(start_distances_y)) {
  min_length <- min(start_distances_y[i], stop_distances_y[i])
  y_values <- c(y_values, min_length)
}

distances_x <- width(exon_ranges)
x_vals <- list()
for (i in 1:length(distances_x)) {
    tmp <- 1:distances_x[i]
    x_vals <- append(x_vals, tmp)
}

