file1 = open('RAVEN/testmodels/Broye/HYMOD/Ost-RAVEN.sh','r')
Lines = file1.readlines()
for line in Lines:
    print(f"f\"{line.strip()}\",")
