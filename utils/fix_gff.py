from optparse import OptionParser

def main():
	parser = OptionParser()
	parser.add_option("-g", "--gff", help="Input GFF (to fix)")
	(args, options) = parser.parse_args()
	f = open(args.gff, 'r')
	ids = []
	parent_ids = []
	for i in f.readlines():
		i = i.rstrip("\n").split("\t")
		att_dict = {}
		att = i[8].split(";")
		if '' in att:
			att.remove('')
		att_order = []
		for a in att:
			k,v = a.split("=")
			att_order.append(k)
			if k in att_dict:
				raise Exception("attribute key seen twice")
			att_dict[k] = v
		if 'ID' in att_dict:
			ids.append(att_dict['ID'])
		if 'Parent' in att_dict:
			#the fix:
			if att_dict['Parent'].endswith('-mRNA-1'):
				att_dict['Parent'] = att_dict['Parent'].replace('-RA-', '-')
			parent_ids.append(att_dict['Parent'])
		#create new attribute line:
		new_att_line = []
		for j in att_order:
			new_att_line.append("%s=%s"%(j, att_dict[j]))
		#print new, fixed GFF line:
		print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], ";". join(new_att_line)))

if __name__ == '__main__':
	main()
