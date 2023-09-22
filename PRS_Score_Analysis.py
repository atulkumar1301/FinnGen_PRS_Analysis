#This code is pipe line for PRS calculation
import subprocess
import os

## Change the directories as per requirement
plink2 = "/Volumes/ATUL_6TB/Tools/./plink2"
plink = "/Volumes/ATUL_6TB/Tools/plink_mac_20221210/./plink"
base_file = "/Volumes/ATUL_6TB/Data/Genetic_Data/BioFINDER_1/GRch38/BioFINDER-1_GWAS_Data_GRCh38"
GWAS_Summary = "/Volumes/ATUL_6TB/Data/GWAS_Summary_Statistics/GWAS_Hearing_Impairment.txt"
pwd = "/Volumes/ATUL_6TB/Work/Projects/Cell_Specific_PRSs/Full_Data/A_beta/"
f_n = "Data.txt" #Clinical information file with demographic and diagnosis for each individual


print ("Extracting GWAS data for Sample")
subprocess.run ([plink2, "--pfile", base_file, "--keep", "Keep_Patient.txt", "--maf", "0.05", "--allow-extra-chr", "--make-bed", "--out", "9_QC_GWAS_data"], cwd = pwd) ## GGenerate a Keep Patient.txt file with FID and IID for all the patients which we need to analyze. Change the maf parameter if needed, but mostly 0.05 is optimal
print ("Extraction of GWAS data for Sample Complete")

print ("Calculating the PCA of genetic Data")
subprocess.run ([plink2, "--bfile", "9_QC_GWAS_data", "--allow-extra-chr", "--pca", "--out", "PCA_FILE"], cwd = pwd)
print ("Calculation of the PCA of genetic Data Complete")

print ("Started the Clumping")
subprocess.run ([plink, "--bfile", "9_QC_GWAS_data", "--allow-extra-chr", "--clump", GWAS_Summary, "--clump-kb", "1000", "--clump-p1", "1", "--clump-p2", "1", "--clump-r2", "0.1", "--out", "Clumped_File"], cwd = pwd)
print ("Finished Clumping")

print ("Calling R for making Pre PRS File")
subprocess.call ("/Volumes/ATUL_6TB/Tools/R_Codes/PRS.R", cwd = pwd)
print ("Generation of Pre PRS File Complete")

### PRS calculation
print ("Started Generating PRS Files")
path = pwd + "PRS/"
C_Main = os.path.exists (path)
if not C_Main:
    os.mkdir (path)
f_m = open (path + "SNP_Extract.txt", 'w', 1)
with open (pwd + "Pre_PRS.txt", 'r') as myFile:
    line = myFile.readline ()
    for line in myFile:
        line_list = line.split("\t")
        f_m.write (str (line_list [2]) + "\n")
SNP_Extract = path + "SNP_Extract.txt"
out_file = path + "10_QC_GWAS_data"
subprocess.run ([plink2, "--bfile", "9_QC_GWAS_data", "--extract", SNP_Extract, "--allow-extra-chr", "--make-pgen", "--out", out_file], cwd = pwd)

l = [0.05, 0.005, 0.0005, 0.00005, 0.000005, 0.0000005, 0.00000005]
for i in l:
    f_m_2 = open(path + "p_value_" + str (i) + ".txt", 'w', 1)
    with open (pwd + "Pre_PRS.txt") as myFile_2:
        line_2 = myFile_2.readline ()
        for line_2 in myFile_2:
            line_list_2 = line_2.split ("\t")
            if float (line_list_2 [5]) <= i:
                f_m_2.write (line_list_2 [2] + "\t" + str(line_list_2 [3]).upper()+ "\t" + line_list_2 [4] + "\n")
for i in l:
    f_i = path + "p_value_" + str(i) + ".txt"
    f_o = path + "p_value_" + str(i)
    subprocess.run ([plink2, "--pfile", "10_QC_GWAS_data", "--score", f_i, "--allow-extra-chr", "--out", f_o], cwd = path)
R_code = subprocess.call (["/Volumes/ATUL_6TB/Tools/R_Codes/PRS_Logistic.R", str(f_n)], cwd = path)
