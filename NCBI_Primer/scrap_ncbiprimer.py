file_name = input("Please enter name for xlsx file to be saved")

import time

from tkinter.filedialog import askopenfilename
filename = askopenfilename()


import selenium
from selenium import webdriver
from selenium.webdriver.common.by import By
driver = webdriver.Safari(executable_path=r'/usr/bin/safaridriver')
driver.get('https://www.ncbi.nlm.nih.gov/tools/primer-blast/')
driver.maximize_window()

import requests
from bs4 import BeautifulSoup

#select upload box
upload_box = driver.find_element_by_name('SEQFILE')
upload_box.send_keys(filename)

get_primer = driver.find_element_by_id('b1')
get_primer.click()

while True:
    page = requests.get(driver.current_url)
    page.content
    soup = BeautifulSoup(page.content, 'html.parser')
    content = soup.get_text()
    if "Graphical view" in content:
        break
    else:
        pass

primer_link = driver.current_url

page = requests.get(primer_link)
page.content

soup = BeautifulSoup(page.content, 'html.parser')

import string
character_list=[]
for each_character in string.ascii_uppercase:
    character_list.append(each_character)

number_list=[]
for a in range (1,11):
    number_list.append(a)
    
from openpyxl import Workbook

workbook = Workbook()
sheet = workbook.active

list_1 = ["Set #", "Primer", "Sequence", "Template strand", "Length", 
"Start", "Stop", "Tm", "GC", "Self_complementary", "Self_3\'_complementary"]
list_2 = ["A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1", "I1", "J1","K1"]

for b in range (1,11):
    cell = "A" + str(2*b)
    sheet[cell] = str(b)
    cell_2 = "B" + str(2*b)
    sheet[cell_2] = "Forward"
    cell_3 = "B" + str(2*b+1)
    sheet[cell_3] = "Reverse"

for number in range (0,11):
    sheet[list_2[number]] = list_1[number]

for i in range (0,10):
    row_cell = []
    forward_cell = []
    reverse_cell = []
    each_row_cell=[]
    for c in range (2,12):
        each_row_cell.append(str(character_list[c]))
    for d in range (0,10):
        forward_cell.append(str(each_row_cell[d])+str(2*(i+1)))
        reverse_cell.append(str(each_row_cell[d])+str(2*(i+1)+1))
    pair_info = soup.find_all(class_="prPairInfo")[i]
    table = pair_info.find('table')
    primer_info = table.find_all('td')
    info_primer = primer_info[0:]
    forward_primer = str(info_primer[0]).strip("<td>").strip("</td>")
    sheet[forward_cell[0]] = str(forward_primer)
    forward_primer_strand = str(info_primer[1]).strip("<td>").strip("</td>")
    sheet[forward_cell[1]] = str(forward_primer_strand)
    forward_primer_length = str(info_primer[2]).strip("<td>").strip("</td>")
    sheet[forward_cell[2]] = str(forward_primer_length)
    forward_primer_start = str(info_primer[3]).strip("<td>").strip("</td>")
    sheet[forward_cell[3]] = str(forward_primer_start)
    forward_primer_stop = str(info_primer[4]).strip("<td>").strip("</td>")
    sheet[forward_cell[4]] = str(forward_primer_stop)
    forward_primer_Tm = str(info_primer[5]).strip("<td>").strip("</td>")
    sheet[forward_cell[5]] = str(forward_primer_Tm)
    forward_primer_GC = str(info_primer[6]).strip("<td>").strip("</td>")
    sheet[forward_cell[6]] = str(forward_primer_GC)
    forward_primer_self_complimentary = str(info_primer[7]).strip("<td>").strip("</td>")
    sheet[forward_cell[7]] = str(forward_primer_self_complimentary)
    forward_primer_self_3_complementary = str(info_primer[8]).strip("<td>").strip("</td>")
    sheet[forward_cell[8]] = str(forward_primer_self_3_complementary)
    reverse_primer = str(info_primer[9]).strip("<td>").strip("</td>")
    sheet[reverse_cell[0]] = str(reverse_primer)
    reverse_primer_strand = str(info_primer[10]).strip("<td>").strip("</td>")
    sheet[reverse_cell[1]] = str(reverse_primer_strand)
    reverse_primer_length = str(info_primer[11]).strip("<td>").strip("</td>")
    sheet[reverse_cell[2]] = str(reverse_primer_length)
    reverse_primer_start = str(info_primer[12]).strip("<td>").strip("</td>")
    sheet[reverse_cell[3]] = str(reverse_primer_start)
    reverse_primer_stop = str(info_primer[13]).strip("<td>").strip("</td>")
    sheet[reverse_cell[4]] = str(reverse_primer_stop)
    reverse_primer_Tm = str(info_primer[14]).strip("<td>").strip("</td>")
    sheet[reverse_cell[5]] = str(reverse_primer_Tm)
    reverse_primer_GC = str(info_primer[15]).strip("<td>").strip("</td>")
    sheet[reverse_cell[6]] = str(reverse_primer_GC)
    reverse_primer_self_complimentary = str(info_primer[16]).strip("<td>").strip("</td>")
    sheet[reverse_cell[7]] = str(reverse_primer_self_complimentary)
    reverse_primer_self_3_complementary = str(info_primer[17]).strip("<td>").strip("</td>")
    sheet[reverse_cell[8]] = str(reverse_primer_self_3_complementary)

sheet.column_dimensions['C'].width = 25

file_name = file_name + ".xlsx"

workbook.save(filename=file_name)
