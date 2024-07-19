###############################################################
# Dockerfile to build container images for SnakePolyester
# Based on naotokubota/polyester:1.0
###############################################################

FROM naotokubota/polyester:1.0

# File Author / Maintainer
LABEL maintainer="Naoto Kubota <naotok@ucr.edu>"

# Install optpaese
RUN R -e "install.packages('optparse', repos='http://cran.rstudio.com/')"

# Launch R
CMD ["R"]
