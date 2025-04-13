# Create a new environment for caching file paths to improve performance
file_path_cache <- new.env(parent = emptyenv())
trait_cache <- new.env(parent = emptyenv())
peaks_cache <- new.env(parent = emptyenv())

