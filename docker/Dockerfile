FROM immcantation/suite:devel
LABEL maintainer="Jason Anthony Vander Heiden [jason.vanderheiden@yale.edu], \
                  Susanna Marquez [susanna.marquez@yale.edu]" \
description="Full suite of the Immcantation framework tools and dependencies."

# Install enchantr
RUN PACKAGE="enchantr" \
    && rm -rf /tmp/${PACKAGE} \
    && git clone https://bitbucket.org/kleinstein/${PACKAGE} /tmp/${PACKAGE} \
    && versions update -n ${PACKAGE} -r /tmp/${PACKAGE} -d immcantation \
    && rinstall -p /tmp/${PACKAGE} \
    && (cd /tmp/${PACKAGE} && builds write -n ${PACKAGE} -v $(git describe --abbrev=12 --always --dirty=+)) \
    && rm -r /tmp/${PACKAGE}

# Set commands
CMD echo -e " Version information:\n"\
            "  versions report\n"\
            "Build stamps:\n"\
            "  builds report\n"\
            "Description of available pipelines:\n"\
            "  pipelines report"
