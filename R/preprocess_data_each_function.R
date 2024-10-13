
#' Preprocess Data for Each Cluster Without Formula
#'
#' This function preprocesses the data from a stars object by subsetting based on a cluster
#' and converting the data to a long format with proper indexing for time and spatial location.
#'
#' @param stars_obj A stars object containing spatial data.
#' @param k The cluster number to subset the data for.
#' @param membership A vector or matrix that defines cluster membership for each region.
#' @param time_var The name of the time variable (default is 'time').
#'
#' @return A long-format data frame with ids for spatial and time indexing, suitable for INLA.
#' @examples
#' preprocess_data_each(stars_obj, k = 1, membership, time_var = "time")
#'
preprocess_data_each <- function(stars_obj, k, membership, time_var = "time") {
  # Check if the input is a stars object
  if (!inherits(stars_obj, "stars")) {
    stop("Input must be a 'stars' object.")
  }

  # Get the indices of the regions that belong to the kth cluster
  cluster_regions <- which(membership == k)

  if (is.vector(membership)) {
    # Subset the stars object to only include the regions in the kth cluster
    stars_obj_subset <- stars_obj[, cluster_regions, , drop = FALSE]
  }

  # Convert the subsetted stars object to a long data frame
  data_long <- as.data.frame(stars_obj_subset, long = TRUE)

  # Calculate the number of time points (nt) and spatial locations (nk) in the cluster
  nt <- length(unique(data_long[[time_var]]))  # number of time points
  nk <- length(cluster_regions)  # number of regions in the kth cluster

  # Create id, idt, and ids ensuring correct indexing based on spatial locations and time
  data_long$id <- seq_len(nt * nk)  # Create unique IDs for each observation
  data_long$idt <- rep(seq_len(nt), each = nk)  # Create time index (idt) for each region
  data_long$ids <- rep(seq_len(nk), time = nt)  # Create spatial index (ids) for each time point

  # Remove the first column (if it contains geometry or spatial information)
  data_inla <- data_long[, -1]

  return(data_inla)  # Return the processed data for the kth cluster
}
