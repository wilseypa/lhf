import sys
import subprocess
import os

proc = subprocess.Popen(["/usr/bin/time", '-v', "./LHFmain/LHF", *sys.argv[1:]], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
stdout, stderr = proc.communicate()
print(stdout.decode())

flag = False
file = os.getcwd() + "/LHFstats.csv"
if not os.path.isfile(file) or os.stat(file).st_size == 0:
	flag = True

with open(file, 'a') as file:
	headers = []
	output = []
	for line in stderr.decode().split('\n'):
		if(line.startswith('\t')):
			line = line.split(": ")
			headers.append(line[0].strip())
			output.append(line[1])
	if flag:
		file.write(','.join(headers) + '\n')
	file.write(','.join(output) + '\n')