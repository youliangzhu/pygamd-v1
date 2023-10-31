import os

# 默认CUDA_PATH环境变量：C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2

# 读取文件列表 begin
def separate_file(filepath):
    tem = []
    for file_name in os.listdir(filepath):
        if file_name.endswith('.dll') and (file_name[:8] == 'cudart64' or file_name[:7] == 'cufft64'):
            # print(file_name)
            input_path = os.path.join(filepath, file_name)
            tem.append(input_path)
    return tem

var = os.getenv('CUDA_PATH')                      # CUDA安装中系统默认设置的环境变量名
cuda_path = var + '\\bin'
poetry_path = os.path.dirname(__file__)           # 获取poetry路径
cuda_path = separate_file(cuda_path)              # cuda文件
poetry_cuda_path = separate_file(poetry_path)     # poetry下的cuda文件

if len(cuda_path) == 0:
    print('未找到CUDA依赖！')
    exit()

if len(poetry_cuda_path) > 0:
    print('pygamd在Win环境下的依赖文件已配置！')
    exit()

for dll in cuda_path:
    os.system(f'copy {dll} {poetry_path}')
print('pygamd在Win环境下依赖文件配置成功！')
exit()