
#--------------------------------------------------------------------------------------------------
# UCSC liftover
#--------------------------------------------------------------------------------------------------

# UCSC liftover chains
mkdir -p data/liftover
wget -P data/liftover http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz  
gunzip data/liftover/*.gz

# download liftOver tool from UCSC:
mkdir -p bin
wget -P bin http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod u+x bin/liftOver


# -------------------------------------------------------------------------------------------------
# Genomic Domains
# -------------------------------------------------------------------------------------------------

mkdir -p data/Dixon2012

# TAD calls of Dixon et al 2012, hESC
wget -P data/Dixon2012 http://chromosome.sdsc.edu/mouse/hi-c/hESC.domain.tar.gz  
# extract and rename
tar xvfz data/Dixon2012/hESC.domain.tar.gz -C data/Dixon2012

cp data/Dixon2012/hESC/combined/total.combined.domain data/Dixon2012/hESC.hg18.bed
rm -r data/Dixon2012/hESC/

# liftover to hg38
bin/liftOver \
        data/Dixon2012/hESC.hg18.bed \
        data/liftover/hg18ToHg38.over.chain \
        data/Dixon2012/hESC.hg18.bed.hg38.bed \
        data/Dixon2012/hESC.hg18.bed.hg38.unmapped.bed







