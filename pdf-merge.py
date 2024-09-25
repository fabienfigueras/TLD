from pypdf import PdfWriter

OutputFile="Abacus-308_GGG_175.pdf"
Str1="Abacus-308_GGG_175_"
#Str2="_2024-09-23.pdf"
Str2=".pdf"
Dist=100
i=0
pdfs = []

while Dist <= 1600:
	FileName= Str1+str(Dist)+Str2
	pdfs.append(FileName)
#	print("FileNAme [",i,"] = ",pdfs[i])
	Dist += 100
	i+=1

# merger = PdfMerger()

print(" Merging pdfs...")

merger = PdfWriter()

for pdf in pdfs:
    merger.append(pdf)

merger.write(OutputFile)
merger.close()
