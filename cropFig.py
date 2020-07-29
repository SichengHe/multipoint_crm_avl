import os
import glob

dir_list = ["./solution/new_case_1/",
    "./solution/new_case_2/",
    "./solution/new_case_3/",
    "./solution/old_case_1/",
    "./solution/old_case_2/",
    "./solution/old_case_3/"]


for i in range(len(dir_list)):
    os.chdir(dir_list[i])
    os.system('ls -l')

    pdfs = glob.glob('./*.pdf')
    for pdf in pdfs:
        os.system('pdfcrop ' + pdf)

    os.chdir("../../")