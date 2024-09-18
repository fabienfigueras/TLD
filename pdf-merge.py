from pypdf import PdfWriter

pdfs = ['Abacus-308_GGG_175_100_2024-09-18.pdf','Abacus-308_GGG_175_200_2024-09-18.pdf','Abacus-308_GGG_175_300_2024-09-18.pdf','Abacus-308_GGG_175_400_2024-09-18.pdf','Abacus-308_GGG_175_500_2024-09-18.pdf','Abacus-308_GGG_175_600_2024-09-18.pdf','Abacus-308_GGG_175_700_2024-09-18.pdf','Abacus-308_GGG_175_800_2024-09-18.pdf','Abacus-308_GGG_175_900_2024-09-18.pdf','Abacus-308_GGG_175_1000_2024-09-18.pdf','Abacus-308_GGG_175_1100_2024-09-18.pdf','Abacus-308_GGG_175_1200_2024-09-18.pdf','Abacus-308_GGG_175_1300_2024-09-18.pdf','Abacus-308_GGG_175_1400_2024-09-18.pdf','Abacus-308_GGG_175_1500_2024-09-18.pdf','Abacus-308_GGG_175_1600_2024-09-18.pdf']

# merger = PdfMerger()

merger = PdfWriter()

for pdf in pdfs:
    merger.append(pdf)

merger.write("Abacus-308_GGG_175.pdf")
merger.close()