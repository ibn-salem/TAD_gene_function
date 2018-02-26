
#===============================================================================
# UCSC liftover
#===============================================================================

# UCSC liftover chains
mkdir -p data/liftover
wget -P data/liftover http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz  
wget -P data/liftover http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz  
gunzip data/liftover/*.gz

# download liftOver tool from UCSC:
mkdir -p bin
wget -P bin http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod u+x bin/liftOver


#===============================================================================
# Genomic Domains
#===============================================================================

#-------------------------------------------------------------------------------
# TADs from Dixen et al. 2012
#-------------------------------------------------------------------------------
mkdir -p data/Dixon2012

# define cell types to download
CELL_TYPES="hESC IMR90"

# iterate over cell types
for CELL in $CELL_TYPES; do
  # TAD calls of Dixon et al 2012, 
  wget -P data/Dixon2012 http://chromosome.sdsc.edu/mouse/hi-c/${CELL}.domain.tar.gz  
  # extract and rename
  tar xvfz data/Dixon2012/${CELL}.domain.tar.gz -C data/Dixon2012
  # rename and remove gunzip archive
  cp data/Dixon2012/${CELL}/combined/total.combined.domain data/Dixon2012/${CELL}.hg18.bed
  rm -r data/Dixon2012/${CELL}/
  # liftover to hg38
  bin/liftOver \
        data/Dixon2012/${CELL}.hg18.bed \
        data/liftover/hg18ToHg38.over.chain \
        data/Dixon2012/${CELL}.hg18.bed.hg38.bed \
        data/Dixon2012/${CELL}.hg18.bed.hg38.unmapped.bed
done

#-------------------------------------------------------------------------------
# TADs from Schmitt et al. 2016
#-------------------------------------------------------------------------------
wget -P data https://ars.els-cdn.com/content/image/1-s2.0-S2211124716314814-mmc4.xlsx

#-------------------------------------------------------------------------------
# TADs from Rao et al. 2014
#-------------------------------------------------------------------------------

mkdir -p data/Rao2014/

# define cell types to download
CELL_TYPES="GM12878_primary+replicate IMR90"

# iterate over cell types
for CELL in $CELL_TYPES; do
  wget -P 'data/Rao2014/' ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_Arrowhead_domainlist.txt.gz
  gunzip 'data/Rao2014/'GSE63525_${CELL}_Arrowhead_domainlist.txt.gz
  # reformat into bed file
  tail -n +2 'data/Rao2014/'GSE63525_${CELL}_Arrowhead_domainlist.txt \
  			|cut -f 1-3 \
  			| sed -e 's/^/chr/' \
  			> 'data/Rao2014/'GSE63525_${CELL}_Arrowhead_domainlist.bed
  # liftover to hg38
  bin/liftOver \
          'data/Rao2014/'GSE63525_${CELL}_Arrowhead_domainlist.bed \
          data/liftover/hg19ToHg38.over.chain \
          'data/Rao2014/'GSE63525_${CELL}_Arrowhead_domainlist.bed.hg38.bed \
          'data/Rao2014/'GSE63525_${CELL}_Arrowhead_domainlist.bed.hg38.unmapped.bed
done





