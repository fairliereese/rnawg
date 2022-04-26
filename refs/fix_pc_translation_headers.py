ifile = open('gencode.v29.pc_translations.fa', 'r')
ofile = open('gencode.v29.pc_translations_short_headers.fa', 'w')
pids = []
curr_pid = ''
for line in ifile:
	if line.startswith('>'):
		line = line.split('|')[0]+'\n'
		curr_pid = line[:-1]
	if curr_pid not in pids:
		ofile.write(line)
	if not line.startswith('>'):
		pids.append(curr_pid)

ifile.close()
ofile.close()
