import os
import shutil
import subprocess


def copy_pht3d_datab(source_folder="caseX"):
    """
    该函数用于将指定案例的 PHREEQC 数据库文件 pht3d_datab.dat 移动到模型文件夹里

    参数:
    - source_folder (str): 源文件夹的路径。
    - destination_folder (str): 目标文件夹的路径，若不存在将自动创建。
    """
    # 定义相对路径
    source_file = os.path.join('.', 'data', source_folder, 'pht3d_datab.dat')
    destination_file = os.path.join('models_folder', source_folder, 'pht3d_datab.dat')
    # 复制并重命名文件
    shutil.copy(source_file, destination_file)
    print(f'The file pht3d_datab.dat has been copied to the model folder.')


def run_pht3d_program(work_dir="caseX", executable="../../bin/pht3dv210.exe", name_file="pht3d.nam"):
    """
    此函数用于在指定的工作目录中运行 PHT3D

    参数:
    - work_dir: 指定PHT3D模拟的工作目录
    - executable: PHT3D程序的执行文件名, 默认为 ../bin/pht3dv210.exe
    - name_file: PHT3D程序的名称文件名, 默认为 pht3d.nam

    异常:
    - FileNotFoundError: 如果工作目录或执行文件不存在，将抛出此异常。
    """
    work_dir = os.path.join("./models_folder", work_dir)
    # D:\Archives\Python\MyPHT3D_model\
    # 检查工作目录是否存在
    if not os.path.exists(work_dir):
        raise FileNotFoundError(f"The specified work directory does not exist: {work_dir}")

    # 检查可执行文件是否存在
    executable_path = os.path.join(work_dir, executable)
    if not os.path.isfile(executable_path):
        raise FileNotFoundError(f"The specified executable file does not exist: {executable_path}")

    # 运行程序并捕获输出
    command = [executable_path, name_file]
    try:
        result = subprocess.run(command, cwd=work_dir, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Command '{e.cmd}' returned non-zero exit status {e.returncode}.")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
    except FileNotFoundError as e:
        print(f"FileNotFoundError: {e}")