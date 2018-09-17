
CPUcount = 6

file = open('Data.txt')
boards = [file.readline().strip().split(' ') for cpu in range(CPUcount)]
print(boards)