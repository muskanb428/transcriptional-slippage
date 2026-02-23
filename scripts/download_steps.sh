#downloading fasta and gtf file

#1. Homo_sapiens
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

#2. Plasmodium_falciparum
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/765/GCF_000002765.6_GCA_000002765/GCF_000002765.6_GCA_000002765_genomic.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/765/GCF_000002765.6_GCA_000002765/GCF_000002765.6_GCA_000002765_genomic.gtf.gz

#3. Arabidopsis_thaliana
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gtf.gz

#4. Trypanosoma_brucei
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/445/GCF_000002445.2_ASM244v1/GCF_000002445.2_ASM244v1_genomic.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/445/GCF_000002445.2_ASM244v1/GCF_000002445.2_ASM244v1_genomic.gtf.gz

#5. Saccharomyces_cerevisiae
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/965/239/625/GCA_965239625.1_Saccharomyces_cerevisiae_NGY10_Assembly/GCA_965239625.1_Saccharomyces_cerevisiae_NGY10_Assembly_genomic.fna.gz

wget https://ftp.ensembl.org/pub/release-115/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gtf.gz

#6. Danio_rerio
wget https://ftp.ensembl.org/pub/release-115/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_rm.primary_assembly.fa.gz

wget https://ftp.ensembl.org/pub/release-115/gtf/danio_rerio/Danio_rerio.GRCz11.115.gtf.gz

#7. Octopus_bimaculoides
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/194/135/GCF_001194135.2_ASM119413v2/GCF_001194135.2_ASM119413v2_genomic.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/194/135/GCF_001194135.2_ASM119413v2/GCF_001194135.2_ASM119413v2_genomic.gtf.gz

#8. Gallus_gallus
 wget https://ftp.ensembl.org/pub/release-115/fasta/gallus_gallus_gca000002315v5/dna/Gallus_gallus_gca000002315v5.GRCg6a.dna.primary_assembly.1.fa.gz
 
wget https://ftp.ensembl.org/pub/release-115/gtf/gallus_gallus_gca000002315v5/Gallus_gallus_gca000002315v5.GRCg6a.115.gtf.gz

############################################################################################################################################################################################################################################
#INDEXING OF FASTA

#1. Homo_sapiens : samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
#2. Arabidopsis_thaliana : samtools faidx GCF_000001735.4_TAIR10.1_genomic.fna 
#3. Danio_rerio : samtools faidx Danio_rerio.GRCz11.dna_rm.primary_assembly.fa
#4. Trypanosoma_brucei : samtools faidx GCF_000002445.2_ASM244v1_genomic.fna
#5. Plasmodium_falciparum : samtools faidx GCF_000002765.6_GCA_000002765_genomic.fna 
#6. Octopus_bimaculoides : samtools faidx GCF_001194135.2_ASM119413v2_genomic.fna 
#7. Saccharomyces_cerevisiae : samtools faidx GCA_965239625.1_Saccharomyces_cerevisiae_NGY10_Assembly_genomic.fna
#8. Gallus_gallus : samtools faidx GCF_000002315.6_GRCg6a_genomic.fna  

