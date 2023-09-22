#! /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/bin/Rscript
##Change the first line as per your machine for the installed location of R and make this code executable using chmod 755
library (data.table)
GWAS_Sum <- fread ("/Volumes/ATUL_6TB/Data/GWAS_Summary_Statistics/GWAS_Hearing_Impairment.txt")
clump <- fread ("Clumped_file.clumped")
write.table (clump [, 3], "clumped_SNP.txt", sep = "\t", quote = FALSE, row.names = FALSE)
Pre_PRS <- merge (clump, GWAS_Sum, by = "SNP")
df_Dup <- Pre_PRS [duplicated(Pre_PRS$SNP)]
df <- Pre_PRS [!duplicated(Pre_PRS$SNP)]
df1 <- subset (df, select = c (CHR_y.y, POS, SNP, EA, Beta, P.y)) # if the Finngen data is in GrCh37
df1 <- subset (df, select = c (CHR_x.y, POS, SNP, EA, Beta, P.y)) # if the Finngen data is in GrCh38
write.table (df1, file = "Pre_PRS.txt", sep = "\t", quote = FALSE, row.names = FALSE)
