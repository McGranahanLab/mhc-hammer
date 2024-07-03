#!/bin/bash -eu

# Declare vars. 
proj_dir=""
bin_dir=""
singularity_image=""
hlahd_downlaod=""

while getopts ":p:h:" opt; do
    case $opt in
        p)
            proj_dir=$OPTARG
            ;;

        h) 
            hlahd_download=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument" >&2
            exit 1
    esac
done
shift $((OPTIND-1))

# check the proj_dir exists
if [ ! -d "$proj_dir" ]; then
    echo "Error: project directory does not exist"
    exit 1
fi

# if hlahd_download is a full path, only keep the filename
if [[ $hlahd_download == */* ]]; then
    hlahd_download=$(basename $hlahd_download)
fi

bin_dir="$proj_dir/bin"
assets_dir="$proj_dir/assets"
hlahd_download="$bin_dir/$hlahd_download"
singularity_image=$proj_dir/singularity_images/mhc_hammer_preprocessing_latest.sif

# check the hlahd_download exists
if [ ! -f "$hlahd_download" ]; then
    echo "Error: hlahd download does not exist"
    exit 1
fi

# check the singularity image exists
if [ ! -f "$singularity_image" ]; then
    echo "Error: singularity image does not exist"
    exit 1
fi

singularity_command="singularity exec -B $proj_dir:$proj_dir $singularity_image"

# move to the bin directory
cd "$bin_dir"

# Install bowtie2
echo "Installing bowtie2 (2.5.1)"
$singularity_command wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.1/bowtie2-2.5.1-linux-x86_64.zip && \
    $singularity_command unzip bowtie2-2.5.1-linux-x86_64.zip && \
    rm bowtie2-2.5.1-linux-x86_64.zip && \
    mv bowtie2-2.5.1-linux-x86_64 bowtie2

# Check installation
if [[ -f bowtie2/bowtie2 && -f bowtie2/bowtie2-build && -f bowtie2/bowtie2-inspect ]]; then
    echo "Files exist. Checking version..."
    bowtie2_version=$(bowtie2/bowtie2 --version 2>&1 | head -n 1 | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    if [[ $bowtie2_version == '2.5.1' ]]; then
        echo "Bowtie2 (2.5.1) installed successfully"
    else
        echo "Unexpected version: $bowtie2_version"
    fi
else
    echo "Bowtie2 installation seems to have failed. Files are missing."
fi

# add bowtie2 to PATH
echo "Adding bowtie2 to PATH"
export PATH="$PATH:$bin_dir/bowtie2"

# now install hlahd
echo "Installing HLA-HD"
$singularity_command tar -xvf "$hlahd_download"
# rm "$hlahd_download"
hlahd_without_extension="${hlahd_download%.tar.gz}"
mv "$hlahd_without_extension" hlahd
cd hlahd

# the following require wget, g++ and bowtie2 to be installed and in the PATH
sh install.sh

# update the IMGT dictionary 
sh ../update.dictionary.alt.sh 