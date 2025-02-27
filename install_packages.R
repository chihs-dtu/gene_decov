# Define a vector with the required package names
required_packages <- c("ggplot2", "dplyr", "tidyr", "RESET", "reticulate", 
                       "pryr")
bio_pacages <- c("fgsea", "AUCell", "doMC", "doRNG", "doSNOW")
# Function to check and install packages
install_packages <- function(packages) {
  # Loop through each package in the list
  for (package in packages) {
    # Check if the package is already installed
    if (!(package %in% installed.packages()[,"Package"])) {
      # If not installed, install it
      message(paste("Installing", package))
      install.packages(package, dependencies = TRUE, 
                      repos = "https://mirrors.dotsrc.org/cran/")
    } else {
      message(paste(package, "is already installed"))
    }
  }
}

# Run the installation function
install_packages(required_packages)
install_packages(bio_pacages)


# Get the list of installed packages and their version numbers
installed_pkg_info <- installed.packages()

# Extract package names and versions
package_versions <- data.frame(
  Package = installed_pkg_info[, "Package"],
  Version = installed_pkg_info[, "Version"],
  stringsAsFactors = FALSE
)


# Write the package names and versions to a text file
write.table(package_versions, file = "installed_packages_versions.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Optionally, print a message
cat("Installed packages and their versions have been saved to 'installed_packages_versions.txt'.\n")


# Optionally, load the installed packages (you can check if they're loaded already)
lapply(required_packages, library, character.only = TRUE)

# # Install python packages
# # Or, use a specific virtualenv
# dirname <- getwd()
# print(dirname)
# virtualenv_create(file.path(dirname,"python_venv"))

# use_virtualenv(file.path(dirname,"python_venv"), required = TRUE)
# py_install("gseapy==1.1.4")
# py_install("rds2py")