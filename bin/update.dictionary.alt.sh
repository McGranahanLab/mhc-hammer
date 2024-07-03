#This shell command updates hla dictionary
#Please ensure the release version at ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/
#Wget program (https://www.gnu.org/software/wget/) is need to run the command.

# if an hla.dat.zip file exists, delete it
if [ -e hla.dat.zip ]; then
rm hla.dat.zip
fi

# if an hla.dat file exists, delete it
if [ -e hla.dat ]; then
rm hla.dat
fi
wget https://github.com/ANHIG/IMGTHLA/raw/3550/hla.dat.zip ## this downloads version 3.55 

unzip -q hla.dat.zip

if [ -e ./bin/create_fasta_from_dat ]; then
:
else
g++ ./src/Create_fasta_from_dat.cpp -O3 -o ./bin/create_fasta_from_dat
fi

if [ -e dictionary ]; then
:
else
mkdir dictionary
fi

mv hla.dat ./dictionary

./bin/create_fasta_from_dat dictionary/hla.dat dictionary/ 150

cd dictionary
sh create_dir.sh
sh move_file.sh
cat full_U_seq.fasta >> all_exon_N150.fasta
cat full_U_seq.convert.list >> convert.list
sh bw_build.sh
cd ../