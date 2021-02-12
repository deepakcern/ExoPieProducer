
def setColumns(df,colmns):
	for col in colmns:
		print "setting coulmn %s to default value 1"%col
		df[col]=1.0

def remove(array,elements):
	for ele in elements:
	    array.remove(ele)
