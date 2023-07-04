import os
import string
import math
import re
import platform
# import logging

curr_path = os.getcwd()
curr_file = os.path.abspath(__file__)

p=curr_file.rfind(".py")
if p!=-1:
	curr_file = curr_file[0:p]

os.system("chmod +x "+curr_file)
os.system("cp "+curr_file+" "+curr_path)
print("copy '"+curr_file+"' to current path: '"+curr_path+"'")
print("load dataTackle successfully!")