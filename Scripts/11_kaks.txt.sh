#----------------------------------------------------------------------------------------
# What: This program combines the kaks files from the same directory and outputs into a new text file
#       The new text file will include the header and will exclude any files that are not in green_anole\tNG format
#       This program will also create csv files which will be used later for running R scripts
#
# Need: axt files made from fasta files
#----------------------------------------------------------------------------------------

awk '!a[$0]++' AAWZ02037698/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/AAWZ02037698.txt
awk '!a[$0]++' AAWZ02040127/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/AAWZ02040127.txt
awk '!a[$0]++' AAWZ02041607/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/AAWZ02041607.txt

# autosomes
awk '!a[$0]++' Chr1/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/Chr1.txt 
awk '!a[$0]++' Chr2/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/Chr2.txt
awk '!a[$0]++' Chr3/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/Chr3.txt 
awk '!a[$0]++' Chr4/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/Chr4.txt
awk '!a[$0]++' Chr5/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/Chr5.txt 
awk '!a[$0]++' Chr6/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/Chr6.txt
awk '!a[$0]++' KaKs_Out/txtfiles/Chr* > KaKs_Out/txtfiles/autosomes.txt

# all scaffolds
awk '!a[$0]++' GL343282/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343282.txt
awk '!a[$0]++' GL343338/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343338.txt
awk '!a[$0]++' GL343364/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343364.txt
awk '!a[$0]++' GL343417/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343417.txt
awk '!a[$0]++' GL343422/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343422.txt
awk '!a[$0]++' GL343423/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343423.txt
awk '!a[$0]++' GL343439/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343439.txt
awk '!a[$0]++' GL343462/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343462.txt
awk '!a[$0]++' GL343516/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343516.txt
awk '!a[$0]++' GL343525/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343525.txt
awk '!a[$0]++' GL343550/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343550.txt
awk '!a[$0]++' GL343588/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343588.txt
awk '!a[$0]++' GL343731/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343731.txt
awk '!a[$0]++' GL343913/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343913.txt
awk '!a[$0]++' GL343947/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL343947.txt
awk '!a[$0]++' GL344042/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL344042.txt
awk '!a[$0]++' GL344393/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL344393.txt
awk '!a[$0]++' GL344496/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL344496.txt
awk '!a[$0]++' GL344539/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/GL344539.txt
awk '!a[$0]++' LGb/8/*.fasta.axt.kaks | grep -v "green_anole	NG" > KaKs_Out/txtfiles/LGb.txt

# scaffolds included in analysis
awk '!a[$0]++' KaKs_Out/txtfiles/LGb.txt KaKs_Out/txtfiles/GL*> KaKs_Out/txtfiles/HypX.txt

# save as csv file

sed 's/	/,/g' KaKs_Out/txtfiles/AAWZ02037698.txt >  KaKs_Out/AAWZ02037698.csv
sed 's/	/,/g' KaKs_Out/txtfiles/AAWZ02040127.txt >  KaKs_Out/AAWZ02040127.csv
sed 's/	/,/g' KaKs_Out/txtfiles/AAWZ02041607.txt >  KaKs_Out/AAWZ02041607.csv
sed 's/	/,/g' KaKs_Out/txtfiles/Chr1.txt >  KaKs_Out/Chr1.csv 
sed 's/	/,/g' KaKs_Out/txtfiles/Chr2.txt >  KaKs_Out/Chr2.csv
sed 's/	/,/g' KaKs_Out/txtfiles/Chr3.txt >  KaKs_Out/Chr3.csv 
sed 's/	/,/g' KaKs_Out/txtfiles/Chr4.txt >  KaKs_Out/Chr4.csv
sed 's/	/,/g' KaKs_Out/txtfiles/Chr5.txt >  KaKs_Out/Chr5.csv 
sed 's/	/,/g' KaKs_Out/txtfiles/Chr6.txt >  KaKs_Out/Chr6.csv
sed 's/	/,/g' KaKs_Out/txtfiles/autosomes.txt > KaKs_Out/autosomes.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343282.txt >  KaKs_Out/GL343282.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343338.txt >  KaKs_Out/GL343338.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343364.txt >  KaKs_Out/GL343364.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343417.txt >  KaKs_Out/GL343417.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343422.txt >  KaKs_Out/GL343422.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343423.txt >  KaKs_Out/GL343423.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343439.txt >  KaKs_Out/GL343439.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343462.txt >  KaKs_Out/GL343462.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343516.txt >  KaKs_Out/GL343516.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343525.txt >  KaKs_Out/GL343525.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343550.txt >  KaKs_Out/GL343550.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343588.txt >  KaKs_Out/GL343588.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343731.txt >  KaKs_Out/GL343731.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343913.txt >  KaKs_Out/GL343913.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL343947.txt >  KaKs_Out/GL343947.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL344042.txt >  KaKs_Out/GL344042.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL344393.txt >  KaKs_Out/GL344393.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL344496.txt >  KaKs_Out/GL344496.csv
sed 's/	/,/g' KaKs_Out/txtfiles/GL344539.txt >  KaKs_Out/GL344539.csv
sed 's/	/,/g' KaKs_Out/txtfiles/LGb.txt >  KaKs_Out/LGb.csv
sed 's/	/,/g' KaKs_Out/txtfiles/HypX.txt > KaKs_Out/HypX.csv

