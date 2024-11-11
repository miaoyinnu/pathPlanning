import matplotlib.pyplot as plt

# 列表中的点坐标
points = [
    (44, 25),  # 点1
    (45, 31),  # 点2
    (48, 38),  # 点3
    (48, 39),  # 点4
    (26, 20),  # 点5
    (27, 18),  # 点6
    (26, 15),  # 点7
    (35, 15),  # 点8
    (42, 13),  # 点9
    (38, 17),  # 点10
    (38, 20),  # 点11
    (42, 23),  # 点12
    (44, 25)   # 回到起点
]

# 分离 X 和 Y 坐标
x_coords = [x for x, y in points]
y_coords = [y for x, y in points]

# 创建绘图
plt.figure(figsize=(10, 8))
plt.plot(x_coords, y_coords, marker='o', linestyle='-', color='b')

# 为每个点添加标注
for i, (x, y) in enumerate(points[:-1]):  # 不包括最后一个重复的起点
    plt.text(x, y, f'{i+1}', fontsize=12, ha='right', va='bottom')

# 设置标题和标签
plt.title('Visual', fontsize=16)
plt.xlabel('X ', fontsize=14)
plt.ylabel('Y ', fontsize=14)
plt.grid(True)

# 设置坐标轴范围
plt.xlim(min(x_coords) - 5, max(x_coords) + 5)
plt.ylim(min(y_coords) - 5, max(y_coords) + 5)

# 显示绘图
plt.show()
