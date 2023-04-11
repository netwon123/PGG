
# 定义需要分类的结构变异类型
categories = ['INS', 'DEL', 'DUP', 'INV', 'BND']

# 打开VCF文件
with open('TT_pop.vcf', 'r') as f:
    # 初始化一个字典来存储结构变异类型和数量
    sv_types = {'unknown': 0}
    # 遍历文件中的每一行
    for line in f:
        # 检查行是否以#开头，如果是，则跳过该行
        if line.startswith('#'):
            continue
        # 将行分割成一个列表
        line_list = line.strip().split('\t')
        # 获取结构变异类型，需要根据vcf修改第一个split的值
        sv_type = line_list[7].split(';')[1].split('=')[1]
        # 将不在分类列表中的结构变异类型归类为unknown
        if sv_type not in categories:
            sv_type = 'unknown'
        # 增加该结构变异类型的数量
        if sv_type in sv_types:
            sv_types[sv_type] += 1
        else:
            sv_types[sv_type] = 1

# 打印每个结构变异类型及其数量
for sv_type, count in sv_types.items():
    print('{}: {}'.format(sv_type, count))
