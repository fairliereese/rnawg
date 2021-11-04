ifile = open('gencode.v29.pc_translations.fa', 'r')
ofile = open('gencode.v29.pc_translations_short_headers.fa', 'w')
for line in ifile:
	if line.startswith('>'):
		line = line.split('|')[0]+'\n'
	ofile.write(line)


ifile.close()
ofile.close()
