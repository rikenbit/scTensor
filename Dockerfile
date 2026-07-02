# Base Image
FROM bioconductor/bioconductor_docker:devel

# Install R Packages
# 'remotes' is required by install_github but is not present in the base image,
# so install it explicitly before use.
RUN R -e "install.packages('remotes'); \
    remotes::install_github('rikenbit/scTensor', \
    upgrade='always', force=TRUE, INSTALL_opts = '--install-tests');\
    tools::testInstalledPackage('scTensor')"
