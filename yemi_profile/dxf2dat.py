file_path = "yemi.dxf"  # Replace with the path to your ASCII file

x = 0
y = 0
z = 0
n = 0

ofile = open("yemi.dat", 'w')

with open(file_path, "r") as file:
    ntot = len(file.readlines())
    print('total number of line in file = %d' % ntot)
file.close()

with open(file_path, "r") as file:
    for line in file:
        n += 1
        obj = line.strip()
        if obj == 'VERTEX':
            while True:
                sline = file.readline().strip()
                if sline == '10':
                    ssline = file.readline().strip()
                    x = float(ssline)
                if sline == '20':
                    ssline = file.readline().strip()
                    y = float(ssline)
                if sline == '30':
                    ssline = file.readline().strip()
                    z = float(ssline)
                    break
            if z > 0:
                ofile.write("%f %f %f \n" % (x, y, z))
                

