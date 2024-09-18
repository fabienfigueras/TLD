#!/usr/bin/python3

import sys
import asyncio
from pyppeteer import launch

# Software introduction for users

h2p_version="1.00"
h2p_year="2024"

print("==================================================")
print(" h2p stands for HTML TO PDF ")
print(" h2p is an Open Source Balistic Software ")
print(" Written in Object Oriented Python3 by Fabien FIGUERAS (he/him)")
print(" v1.00 was released in 2024 ")
print(" Current Version is v",h2p_version,h2p_year)
print(" Call example h2p ./h2p-vxyz.py to get this message")
print("")
print(" Where param1 is the name of the HTML file to be converted into pdf")
print(" Where param2 is the name of the pdf file ")
print("")
print(" Sources available in GitHub : https://github.com/fabienfigueras/TLD")
print("==================================================")

FName_H=sys.argv[1]

print("HTML File Name :",FName_H)

FName_P=sys.argv[2]

print("PDF File Name :",FName_P)

async def generate_pdf_from_html(html_content, pdf_path):
    browser = await launch()
    page = await browser.newPage()

    await page.setContent(html_content)

    await page.pdf({'path': pdf_path, 'format': 'A4'})

    await browser.close()


F_H=open(FName_H,'r')
# HTML content
html_content = F_H.read()

# print(html_content)

asyncio.get_event_loop().run_until_complete(generate_pdf_from_html(html_content, FName_P))

F_H.close()