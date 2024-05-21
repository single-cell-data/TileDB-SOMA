#!/bin/bash
# -*- sh-basic-offset: 4; sh-indentation: 4 -*-
#
# Bootstrap a CI environment for R
#
# Originally written by Craig Citro and contributors, see https://github.com/craigcitro/r-travis
# Later forked and maintained by Dirk Eddelbuettel, see https://github.com/eddelbuettel/r-ci
#
# Released under Apache 2.0 License

# Contains some extraneous code which could be removed but has supported (tens of) thousands of CI runs

set -e
# Comment out this line for quieter output:
# Or ratherm set it for a lot noisier output
# set -x

# Where do we run? OS can be 'Linux' or 'Darwin', ARCH can be 'x86_64' or 'arm64'
OS=$(uname -s)
ARCH=$(uname -m)

# Default CRAN repo (use the CDN) and R verssion
CRAN=${CRAN:-"https://cloud.r-project.org"}
RVER=${RVER:-"4.3.3"}

## Optional drat repos, unset by default
DRAT_REPOS=${DRAT_REPOS:-""}

## Optional BSPM use, defaults to true for r2u
USE_BSPM=${USE_BSPM:-"TRUE"}

## Optional additional PPAs, unset by default
ADDED_PPAS=${ADDED_PPAS:-""}

## Optional trimming of extra apt source list entry, defaults to false
TRIM_APT_SOURCES=${TRIM_APT_SOURCES:-"TRUE"}

## Optional setting of type argument in covr::coverage() call below, defaults to "tests"
COVERAGE_TYPE=${COVERAGE_TYPE:-"tests"}
COVERAGE_FLAGS=${COVERAGE_FLAGS:-"r"}
COVERAGE_TOKEN=${COVERAGE_TOKEN:-""}
COVERAGE_PATH=${COVERAGE_PATH:-"apis/r"}

## Optional 'catchsegv' frontend on Linux
CATCHSEGV=${CATCHSEGV:-"FALSE"}

R_BUILD_ARGS=${R_BUILD_ARGS-"--no-build-vignettes --no-manual"}
R_CHECK_ARGS=${R_CHECK_ARGS-"--no-manual --as-cran"}
R_CHECK_INSTALL_ARGS=${R_CHECK_INSTALL_ARGS-"--install-args=--install-tests"}
_R_CHECK_TESTS_NLINES_=0

ShowBanner() {
    echo ""
    echo "r-ci: Portable CI for R at Travis, GitHub Actions, Azure, ..."
    echo ""
    echo "Current Ubuntu distribution per 'lsb_release': '$(lsb_release -ds)' aka '$(lsb_release -cs)'"
    echo ""
}

Bootstrap() {
    SetRepos

    if [[ "Darwin" == "${OS}" ]]; then
        BootstrapMac
    elif [[ "Linux" == "${OS}" ]]; then
        BootstrapLinux
    else
        echo "Unknown OS: ${OS}"
        exit 1
    fi

    if ! (test -e .Rbuildignore && grep -q 'travis-tool' .Rbuildignore); then
        echo '^travis-tool\.sh$' >>.Rbuildignore
    fi
    if ! (test -e .Rbuildignore && grep -q 'run.sh' .Rbuildignore); then
        echo '^run\.sh$' >> .Rbuildignore
    fi
    if ! (test -e .Rbuildignore && grep -q 'travis_wait' .Rbuildignore); then
        echo '^travis_wait_.*\.log$' >> .Rbuildignore
    fi

    # Make sure unit test package (among testthat, tinytest, RUnit) installed
    EnsureUnittestRunner

    # Report version
    Rscript -e 'sessionInfo()'
}

SetRepos() {
    echo "local({" >> ~/.Rprofile
    echo "   r <- getOption(\"repos\");" >> ~/.Rprofile
    echo "   r[\"CRAN\"] <- \"${CRAN}\"" >> ~/.Rprofile
    for d in ${DRAT_REPOS}; do
        echo "   r[\"${d}\"] <- \"https://${d}.github.io/drat\"" >> ~/.Rprofile
    done
    echo "   options(repos=r)" >> ~/.Rprofile
    echo "})" >> ~/.Rprofile
}

InstallPandoc() {
    ## deprecated 2020-Sep
    echo "Deprecated"
}

BootstrapLinux() {
    ## Check for sudo_release and install if needed
    test -x /usr/bin/sudo || apt-get install -y --no-install-recommends sudo
    ## Hotfix for key issue
    echo 'Acquire::AllowInsecureRepositories "true";' | sudo tee /etc/apt/apt.conf.d/90local-secure >/dev/null

    ## Check for lsb_release and install if needed
    test -x /usr/bin/lsb_release || sudo apt-get install -y --no-install-recommends lsb-release
    ## Check for add-apt-repository and install if needed, using a fudge around the (manual) tz config dialog
    test -x /usr/bin/add-apt-repository || \
        (echo 12 > /tmp/input.txt; echo 5 >> /tmp/input.txt; sudo apt-get install -y tzdata < /tmp/input.txt; sudo apt-get install -y --no-install-recommends software-properties-common)

    ShowBanner

    ## If opted in, trim apt sources
    if [[ "${TRIM_APT_SOURCES}" != "FALSE" ]]; then
        sudo rm -vf /etc/apt/sources.list.d/*.list
    fi

    ## from r2u setup script
    sudo apt update -qq && sudo apt install --yes --no-install-recommends wget ca-certificates dirmngr gnupg gpg-agent
    wget -q -O- https://eddelbuettel.github.io/r2u/assets/dirk_eddelbuettel_key.asc | sudo tee -a /etc/apt/trusted.gpg.d/cranapt_key.asc
    echo "deb [arch=amd64] https://r2u.stat.illinois.edu/ubuntu $(lsb_release -cs) main" | sudo tee -a /etc/apt/sources.list.d/cranapt.list
    wget -q -O- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc  | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
    echo "deb [arch=amd64] https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" | sudo tee -a /etc/apt/sources.list.d/cran_r.list
    echo "Package: *" | sudo tee -a /etc/apt/preferences.d/99cranapt
    echo "Pin: release o=CRAN-Apt Project" | sudo tee -a /etc/apt/preferences.d/99cranapt
    echo "Pin: release l=CRAN-Apt Packages" | sudo tee -a /etc/apt/preferences.d/99cranapt
    echo "Pin-Priority: 700" | sudo tee -a /etc/apt/preferences.d/99cranapt


    ## Set up our CRAN mirror.
    ## Get the key if it is missing
    if ! test -f /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc; then
       wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
    fi
    ## Add the repo
    ## need pinning to ensure repo sorts higher, note we also pin r2u
    #echo "Package: *" | sudo tee /etc/apt/preferences.d/c2d4u-pin >/dev/null
    #echo "Pin: release o=LP-PPA-c2d4u.team-c2d4u4.0+" | sudo tee -a /etc/apt/preferences.d/c2d4u-pin >/dev/null
    #echo "Pin-Priority: 600" | sudo tee -a /etc/apt/preferences.d/c2d4u-pin >/dev/null
    ## now add repo (and update index)
    # already done above:  sudo add-apt-repository -y "deb ${CRAN}/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

    # Add marutter's c2d4u repository.
    # R 4.0 (not needed as CRAN current) and c2d4u/4.0 variant as backup
    #sudo add-apt-repository -y "ppa:marutter/rrutter4.0"
    #sudo add-apt-repository -y "ppa:c2d4u.team/c2d4u4.0+"

    ## Added PPAs, if given
    if [[ "${ADDED_PPAS}" != "" ]]; then
        for ppa in "${ADDED_PPAS}"; do
            sudo add-apt-repository -y "${ppa}"
        done
    fi


    # Update after adding all repositories.  Retry several times to work around
    # flaky connection to Launchpad PPAs.
    Retry sudo apt-get update -qq

    # Install an R development environment. qpdf is also needed for
    # --as-cran checks:
    #   https://stat.ethz.ch/pipermail/r-help//2012-September/335676.html
    # May 2020: we also need devscripts for checkbashism
    # Sep 2020: add bspm, remotes
    PKGS="r-base-dev r-recommended qpdf devscripts r-cran-bspm r-cran-remotes"

    # If 'catchsegv' requested intall glibc-tools
    if [[ "${CATCHSEGV}" != "FALSE" ]]; then
        PKGS="${PKGS} glibc-tools"
    fi

    Retry sudo apt-get install -y --no-install-recommends ${PKGS}

    #sudo cp -ax /usr/lib/R/site-library/littler/examples/{build.r,check.r,install*.r,update.r} /usr/local/bin
    ## for now also from littler from GH
    #sudo install.r remotes
    #sudo installGithub.r eddelbuettel/littler
    #sudo cp -ax /usr/local/lib/R/site-library/littler/examples/{check.r,install*.r} /usr/local/bin

    # Default to no recommends
    echo 'APT::Install-Recommends "false";' | sudo tee /etc/apt/apt.conf.d/90local-no-recommends >/dev/null

    # Change permissions for /usr/local/lib/R/site-library
    # This should really be via 'staff adduser travis staff'
    # but that may affect only the next shell
    sudo chmod 2777 /usr/local/lib/R /usr/local/lib/R/site-library

    # Process options
    BootstrapLinuxOptions
}

BootstrapLinuxOptions() {
    if [[ -n "$BOOTSTRAP_LATEX" ]]; then
        # We add a backports PPA for more recent TeX packages.
        # sudo add-apt-repository -y "ppa:texlive-backports/ppa"
        Retry sudo apt-get install -y --no-install-recommends \
            texlive-base texlive-latex-base \
            texlive-fonts-recommended texlive-fonts-extra \
            texlive-extra-utils texlive-latex-recommended texlive-latex-extra \
            texinfo lmodern
        # no longer exists: texlive-generic-recommended
    fi
    #if [[ -n "$BOOTSTRAP_PANDOC" ]]; then
    #    InstallPandoc 'linux/debian/x86_64'
    #fi
    if [[ "${USE_BSPM}" != "FALSE" ]]; then
        echo "options(bspm.sudo = TRUE)" | sudo tee --append /etc/R/Rprofile.site >/dev/null
        echo "suppressMessages(bspm::enable())" | sudo tee --append /etc/R/Rprofile.site >/dev/null
        echo "options(bspm.version.check=FALSE)" | sudo tee --append /etc/R/Rprofile.site >/dev/null
    fi
}

BootstrapMac() {
    # Install from latest CRAN binary build for OS X (given ${ARCH} from 'uname -m')
    wget ${CRAN}/bin/macosx/big-sur-${ARCH}/base/R-${RVER}-${ARCH}.pkg  -O /tmp/R-latest.pkg

    echo "Installing macOS binary package for R on ${ARCH}"
    sudo installer -pkg "/tmp/R-latest.pkg" -target /
    rm "/tmp/R-latest.pkg"

    # Process options
    BootstrapMacOptions

    # Default packages
    sudo Rscript -e 'install.packages(c("remotes"))'
}

BootstrapMacOptions() {
    if [[ -n "$BOOTSTRAP_LATEX" ]]; then
        # TODO: Install MacTeX.pkg once there's enough disk space
        MACTEX=BasicTeX.pkg
        wget http://ctan.math.utah.edu/ctan/tex-archive/systems/mac/mactex/$MACTEX -O "/tmp/$MACTEX"

        echo "Installing OS X binary package for MacTeX"
        sudo installer -pkg "/tmp/$MACTEX" -target /
        rm "/tmp/$MACTEX"
        # We need a few more packages than the basic package provides; this
        # post saved me so much pain:
        #   https://stat.ethz.ch/pipermail/r-sig-mac/2010-May/007399.html
        sudo tlmgr update --self
        sudo tlmgr install inconsolata upquote courier courier-scaled helvetic
    fi
    #if [[ -n "$BOOTSTRAP_PANDOC" ]]; then
    #    InstallPandoc 'mac'
    #fi
}

EnsureDevtools() {
    ## deprecated 2020-Sep
    echo "Deprecated"
    #if ! Rscript -e 'if (!("devtools" %in% rownames(installed.packages()))) q(status=1)' ; then
    #    # Install devtools and testthat.
    #    RBinaryInstall devtools testthat
    #fi
}

EnsureUnittestRunner() {
    sudo Rscript -e 'dcf <- read.dcf(file="DESCRIPTION")[1,]; if ("Suggests" %in% names(dcf)) { sug <- dcf[["Suggests"]]; pkg <- do.call(c, sapply(c("testthat", "tinytest", "RUnit"), function(p, sug) if (grepl(p, sug)) p else NULL, sug, USE.NAMES=FALSE)); if (!is.null(pkg)) install.packages(pkg) }'
}

InstallIfNotYetInstalled() {
    res=$(Rscript -e 'if (requireNamespace(commandArgs(TRUE), quietly=TRUE)) cat("YES") else cat("NO")' "$1")
    if [[ "${res}" != "YES" ]]; then
        sudo Rscript -e 'install.packages(commandArgs(TRUE))' "$1"
    fi
}

AptGetInstall() {
    if [[ "Linux" != "${OS}" ]]; then
        echo "Wrong OS: ${OS}"
        exit 1
    fi

    if [[ "" == "$*" ]]; then
        echo "No arguments to aptget_install"
        exit 1
    fi

    echo "Installing apt package(s) $@"
    Retry sudo apt-get -y --no-install-recommends --allow-unauthenticated install "$@"
}

DpkgCurlInstall() {
    if [[ "Linux" != "${OS}" ]]; then
        echo "Wrong OS: ${OS}"
        exit 1
    fi

    if [[ "" == "$*" ]]; then
        echo "No arguments to dpkgcurl_install"
        exit 1
    fi

    echo "Installing remote package(s) $@"
    for rf in "$@"; do
        curl -OL ${rf}
        f=$(basename ${rf})
        sudo dpkg -i ${f}
        rm -v ${f}
    done
}

RInstall() {
    if [[ "" == "$*" ]]; then
        echo "No arguments to r_install"
        exit 1
    fi

    echo "Installing R package(s): $@"
    sudo Rscript -e 'install.packages(commandArgs(TRUE))' "$@"
}

BiocInstall() {
    if [[ "" == "$*" ]]; then
        echo "No arguments to bioc_install"
        exit 1
    fi

    echo "Installing R Bioconductor package(s): $@"
    Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")'
    Rscript -e "BiocManager::install(commandArgs(TRUE), update=FALSE)" "$@"
}

RBinaryInstall() {
    if [[ -z "$#" ]]; then
        echo "No arguments to r_binary_install"
        exit 1
    fi

    if [[ "Linux" != "${OS}" ]] || [[ -n "${FORCE_SOURCE_INSTALL}" ]]; then
        echo "Fallback: Installing from source"
        RInstall "$@"
        return
    fi

    echo "Installing *binary* R packages: $*"
    r_packages=$(echo $* | tr '[:upper:]' '[:lower:]')
    r_debs=$(for r_package in ${r_packages}; do echo -n "r-cran-${r_package} "; done)

    AptGetInstall ${r_debs}
}

InstallGithub() {
    #EnsureDevtools
    sudo Rscript -e 'remotes::install_github(commandArgs(TRUE))' "$@"
}

InstallDeps() {
    #EnsureDevtools
    sudo Rscript -e 'remotes::install_deps(".")'
}

InstallDepsAndSuggests() {
    sudo Rscript -e 'options(timeout = max(300, getOption("timeout"))); remotes::install_deps(".", dependencies=TRUE)'
}

DumpSysinfo() {
    echo "Dumping system information."
    Rscript -e '.libPaths(); sessionInfo(); installed.packages()'
}

DumpLogsByExtension() {
    if [[ -z "$1" ]]; then
        echo "dump_logs_by_extension requires exactly one argument, got: $@"
        exit 1
    fi
    extension=$1
    shift

    # Note: this script may be run from the apis/r subdirectory or from the
    # package base directory. This runs in CI from a clean image so in the
    # (hypothetical) case of multiple log files, we err on the side of providing
    # the operator with more inforation rather than less.
    found="false"
    for package in $(find . -name "*.Rcheck" -type d); do
        for name in $(find "${package}" -type f -name "*${extension}"); do
            found="true"
            echo "::group::Package ${package} filename ${name}" # GHA magic
            cat ${name}
            echo "::endgroup::"
        done
    done
    if [[ "$found" = "false" ]]; then
        echo "Could not find any package Rcheck directories; skipping log dump."
        exit 0
    fi
}

DumpLogs() {
    echo "Dumping test execution logs."
    DumpLogsByExtension "out"
    DumpLogsByExtension "log"
    DumpLogsByExtension "fail"
}

Coverage() {
    echo "Running Code Coverage analysis via the covr package"

    ## assumes that the Rutter PPAs are in fact known, which is a given here
    AptGetInstall r-cran-covr

    # COVERAGE_TYPE=${COVERAGE_TYPE:-"tests"}
    # COVERAGE_FLAGS=${COVERAGE_FLAGS:-"r"}
    # COVERAGE_TOKEN=${COVERAGE_TOKEN:-""}
    # COVERAGE_PATH=${COVERAGE_PATH:-"apis/r"}
    if [[ "${CATCHSEGV}" != "FALSE" ]] && [[ "Linux" == "${OS}" ]]; then
        COVR="true" catchsegv Rscript -e 'setwd("apis/r"); library(covr); res <- package_coverage(relative_path="../..", line_exclusion=list("apis/r/src/nanoarrow.c", "apis/r/src/nanoarrow.h", "apis/r/R/roxygen.R"), quiet=FALSE); print(res); codecov(coverage=res, token="", flags="r")'
    else
        COVR="true" Rscript -e 'setwd("apis/r"); library(covr); res <- package_coverage(relative_path="../..", line_exclusion=list("apis/r/src/nanoarrow.c", "apis/r/src/nanoarrow.h", "apis/r/R/roxygen.R"), quiet=FALSE); print(res); codecov(coverage=res, token="", flags="r")'
    fi
}

RunTests() {
    echo "Building with: R CMD build ${R_BUILD_ARGS}"
    R CMD build ${R_BUILD_ARGS} .
    # We want to grab the version we just built.
    FILE=$(ls -1t *.tar.gz | head -n 1)

    echo "Testing with: R CMD check \"${FILE}\" ${R_CHECK_ARGS} ${R_CHECK_INSTALL_ARGS}"
    _R_CHECK_CRAN_INCOMING_=${_R_CHECK_CRAN_INCOMING_:-FALSE}
    if [[ "$_R_CHECK_CRAN_INCOMING_" == "FALSE" ]]; then
        echo "(CRAN incoming checks are off)"
    fi
    if [[ "${CATCHSEGV}" != "FALSE" ]] && [[ "Linux" == "${OS}" ]]; then
        _R_CHECK_CRAN_INCOMING_=${_R_CHECK_CRAN_INCOMING_} catchsegv R CMD check "${FILE}" ${R_CHECK_ARGS} ${R_CHECK_INSTALL_ARGS}
    else
        _R_CHECK_CRAN_INCOMING_=${_R_CHECK_CRAN_INCOMING_} R CMD check "${FILE}" ${R_CHECK_ARGS} ${R_CHECK_INSTALL_ARGS}
    fi

    if [[ -n "${WARNINGS_ARE_ERRORS}" ]]; then
        if DumpLogsByExtension "00check.log" | grep -q WARNING; then
            echo "Found warnings, treated as errors."
            echo "Clear or unset the WARNINGS_ARE_ERRORS environment variable to ignore warnings."
            exit 1
        fi
    fi
}

Retry() {
    if "$@"; then
        return 0
    fi
    for wait_time in 5 20 30 60; do
        echo "Command failed, retrying in ${wait_time} ..."
        sleep ${wait_time}
        if "$@"; then
            return 0
        fi
    done
    echo "Failed all retries!"
    exit 1
}

ShowHelpAndExit() {
    echo "Usage: run.sh COMMAND"
    echo "Derived from the venerable r-travis project, and still maintained lovingly by @eddelbuettel."
    echo "See https://eddelbuettel.github.io/r-ci for more."
    exit 0
}

COMMAND=$1
#echo "Running command: ${COMMAND}"
shift
case $COMMAND in
    ##
    ## Bootstrap a new core system
    "bootstrap")
        Bootstrap
        ;;
    ## Code coverage via covr.io
    "coverage")
        Coverage
        ;;
    ##
    ## Ensure devtools is loaded (implicitly called)
    "install_devtools"|"devtools_install")
        EnsureDevtools
        ;;
    ##
    ## Install a binary deb package via apt-get
    "install_aptget"|"aptget_install")
        AptGetInstall "$@"
        ;;
    ##
    ## Install a binary deb package via a curl call and local dpkg -i
    "install_dpkgcurl"|"dpkgcurl_install")
        DpkgCurlInstall "$@"
        ;;
    ##
    ## Install an R dependency from CRAN
    "install_r"|"r_install")
        RInstall "$@"
        ;;
    ##
    ## Install an R dependency from Bioconductor
    "install_bioc"|"bioc_install")
        BiocInstall "$@"
        ;;
    ##
    ## Install an R dependency as a binary (via c2d4u PPA)
    "install_r_binary"|"r_binary_install")
        RBinaryInstall "$@"
        ;;
    ##
    ## Install a package from github sources
    "install_github"|"github_package")
        InstallGithub "$@"
        ;;
    ##
    ## Install package dependencies from CRAN
    "install_deps")
        InstallDeps
        ;;
    ##
    ## Install package dependencies and suggests from CRAN
    "install_all")
        InstallDepsAndSuggests
        ;;
    ##
    ## Run the actual tests, ie R CMD check
    "run_tests")
        RunTests
        ;;
    ##
    ## Dump information about installed packages
    "dump_sysinfo")
        DumpSysinfo
        ;;
    ##
    ## Dump build or check logs
    "dump_logs")
        DumpLogs
        ;;
    ##
    ## Dump selected build or check logs
    "dump_logs_by_extension")
        DumpLogsByExtension "$@"
        ;;
    ##
    ## Help
    "help")
        ShowHelpAndExit
        ;;
esac
