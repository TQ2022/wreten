#!/bin/bash

# 获取当前目录
currentDir=$(pwd)

# 遍历当前目录下的所有文件和目录
for file in $currentDir/*; do
    # 检查文件扩展名是否为xmf、h5或res
    if [[ $file == *.xmf || $file == *.h5 || $file == *.res ]]; then
        # 删除文件
        rm "$file"
        echo "已删除文件: $file"
    fi
done