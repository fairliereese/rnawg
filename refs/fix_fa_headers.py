f = 'hg38_sirv4_ercc.fa'

ifile = open(f, 'r')
ofile = open('temp', 'w')

for line in ifile:
    if line.startswith('>'):
        line = line.split(' ')[0]
    ofile.write(line)

ifile.close()
ofile.close()
