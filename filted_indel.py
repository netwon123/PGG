import sys
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('#'):
            print(line.strip())
            continue
        fields = line.strip().split('\t')

        if abs(len(fields[3]) - len(fields[4])) > 49:
            continue
        print(line.strip())
