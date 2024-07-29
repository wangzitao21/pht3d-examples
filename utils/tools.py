import os
# import glob
import shutil
import numpy as np

def clear_directory(directory_path="models_folder"):
    """
    删除指定文件夹里的所有子文件夹和文件，但保留.gitignore文件。

    参数:
    directory_path (str): 目标文件夹的路径。
    """
    if not os.path.exists(directory_path):
        print(f"路径 {directory_path} 不存在。")
        return

    # 遍历目录中的所有文件和子文件夹
    for filename in os.listdir(directory_path):
        # 如果文件名是.gitignore，则跳过
        if filename == ".gitignore":
            continue

        file_path = os.path.join(directory_path, filename)
        try:
            # 如果是文件，则删除
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            # 如果是目录，则删除目录及其内容
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"删除文件或文件夹 {file_path} 时发生错误: {e}")


def calculate_MAS_sums(filename, if_print=False):
    """
    计算并返回MAS文件中每列数据的总和。
    
    参数:
    filename: 字符串, 表示MAS文件的路径。
    if_print: 布尔值, 决定是否打印计算结果。
    
    返回:
    一个字典, 包含每列数据的总和。
    """
    # 定义文件名
    filename = filename

    # 初始化列名和数据
    columns = ["TIME (D)", "TOTAL IN (KG)", "TOTAL OUT (KG)", "SOURCES (KG)", "SINKS (KG)", "NET MASS FROM FLUID-STORAGE (KG)", "TOTAL MASS IN AQUIFER (KG)", "DISCREPANCY (TOTAL IN-OUT)", "DISCREPANCY (ALTERNATIVE)"]
    data = {col: [] for col in columns}

    # 读取文件
    with open(filename, 'r') as file:
        lines = file.readlines()

    # 跳过前两行（表头和单位）
    for line in lines[2:]:
        values = line.split()
        for i, col in enumerate(columns):
            data[col].append(float(values[i]))

    # 计算每列的和
    sums = {col: np.sum(data[col]) for col in columns}

    if if_print:
        print("Sum of each column:")
        for col, sum_value in sums.items():
            print(f"{col}: {sum_value:.2f}")

    return sums