import pandas as pd
import matplotlib.pyplot as plt

# 读取RES文件
data = pd.read_csv('water_retention.res', sep='\\s+')

# 提取Pc和Sr列的数据
Pc = data['Pc']
Sr = data['Sr']

# 设定N的值
N =0  # 设定要去掉的前N个点

# 过滤数据
filtered_Pc = Pc[N:]
filtered_Sr = Sr[N:]

# 创建画布和子图
fig, ax = plt.subplots(figsize=(8, 6))

# 绘制散点图和曲线
ax.scatter(filtered_Sr[::5], filtered_Pc[::5], c='b', label='Sampled Data')
ax.plot(filtered_Sr, filtered_Pc, c='r', label='Curve')

# 设置标题和标签
ax.set_title('Sr vs. Pc', fontsize=16)
ax.set_xlabel('Sr', fontsize=12)
ax.set_ylabel('Pc', fontsize=12)

# 添加网格线和图例
ax.grid(True, linestyle='--', alpha=0.5)
ax.legend()

# 调整坐标轴刻度字体大小
ax.tick_params(axis='both', labelsize=10)

# 调整图像边距
plt.tight_layout()

# 保存图像为PNG文件
plt.savefig('plot.png', dpi=300, format='png')

# 显示图形
plt.show()
