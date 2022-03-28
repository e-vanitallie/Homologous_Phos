#!/usr/bin/env python3

import pandas as pd

def crossref_human_database(myfile)

    df = pd.read_csv(myfile,sep='\t',skiprows=(0,1,2),header=(0))

    print(df.head())
