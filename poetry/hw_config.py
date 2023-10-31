import os

poetry_path = os.path.dirname(__file__)           # 获取poetry路径
print(poetry_path)

if os.system('rm -rf '+poetry_path+'/cu_gala.so'):
    print("Failed for deleting cu_gala library!!")
if not os.system('cp '+poetry_path+'/hw_cu_gala.so '+poetry_path+'/cu_gala.so'):
    print("Successful for adding HW cu_gala library!!")
else:
    print("Failed for adding HW cu_gala library!!")
    
    
if os.system('rm -rf '+poetry_path+'/molgen.so'):
    print("Failed for deleting molgen library!!")
if not os.system('cp '+poetry_path+'/hw_molgen.so '+poetry_path+'/molgen.so'):
    print("Successful for adding HW molgen library!!")
else:
    print("Failed for adding HW molgen library!!")
    
    
if os.system('rm -rf '+poetry_path+'/dataTackle'):
    print("Failed for deleting dataTackle!!")
if not os.system('cp '+poetry_path+'/hw_dataTackle '+poetry_path+'/dataTackle'):
    print("Successful for adding HW dataTackle!!")
else:
    print("Failed for adding HW dataTackle!!")