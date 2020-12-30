import csv
import xml.etree.ElementTree as ET

#element tree object
tree = ET.parse("./catplus.20201101.xml")

#root element
root = tree.getroot()

#loop over each record to find information
title = []
abbre = []
for record in root.findall('NLMCatalogRecord'):
    title_node = record.find('TitleMain')
    full_title = title_node.find("Title").text
    abbreviation_node = record.find('MedlineTA')
    if abbreviation_node is not None:
        abbreviation = abbreviation_node.text
        title.append(full_title)
        abbre.append(abbreviation)

with open('./journal_abbre.txt',"w") as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(title,abbre))




