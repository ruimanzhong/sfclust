
library(testthat)
library(stars)

# Test preprocess_data_each function without formula
test_that('preprocess_data_each works for vector membership without formula', {
  coords <- st_sfc(st_point(c(0, 0)), st_point(c(1, 1)), st_point(c(2, 2)))

  # Time variable
  time_var <- as.POSIXct(c('2022-01-01', '2022-01-02'))

  # Case, temperature, precipitation, population, and expected_cases data arrays
  case <- array(runif(6, 50, 100), dim = c(3, 2))  # 3 spatial points x 2 time points
  temperature <- array(runif(6, 10, 20), dim = c(3, 2))
  precipitation <- array(runif(6, 0, 50), dim = c(3, 2))

  # Combine arrays into a stars object
  stars_obj <- st_as_stars(case = case, temperature = temperature, precipitation = precipitation)

  # Assign dimensions (spatial and time)
  st_dimensions(stars_obj) <- st_dimensions(geometry =coords, time = time_var)


  membership <- c(1, 1, 2)

  # Test the function with vector membership
  processed_data <- preprocess_data_each(stars_obj, k = 1, membership)

  # 1. Check the number of rows
  expect_equal(nrow(processed_data), 4)  # 2 time points * 2 regions in cluster 1

  # 2. Check the ids (spatial) and idt (time) indexes
  expect_equal(processed_data$idt, rep(1:2, each = 2))  # Time index (2 time points, repeated for each region)
  expect_equal(processed_data$ids, rep(1:2, 2))  # Spatial index (2 regions in the cluster)


})
