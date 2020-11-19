import pandas as pd
import sys
#Bigger dataframe Smakker-dataframes' name

test = pd.read_csv(str(sys,argv[1]))
test[(test['Names'].isin(to_iterate))]

test.to_csv(str(sys.argv[2])+"_all.csv", row.names = False)