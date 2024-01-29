import itertools

# 读取文件名列表
with open('path.txt', 'r') as file:
    filenames = file.read().splitlines()

# 生成两两配对的组合
pairs = list(itertools.combinations(filenames, 2))

# 将结果写入txt文件
with open('pairs.txt', 'w') as file:
    for pair in pairs:
        file.write(pair[0] + ' ' + pair[1] + '\n')
        file.write(pair[1] + ' ' + pair[0] + '\n')
