#!/usr/bin/python

from collections import defaultdict
from pyfaidx import Fasta
import sys
import argparse, re
import ctypes
import struct
from optparse import OptionParser
import csv
import pandas as pd
import glob
import numpy as np

import pandas as pd

def create_aspergillus_pooled_category(df):
	aspergillus_tests = [
    	"Aspergillus_flavus",
    	"Aspergillus_fumigatus",
    	"Aspergillus_niger",
    	"Aspergillus_terreus",
    	"Aspergillus_spp"
	]

	# Filter and copy
	df_aspergillus = df[df["mgcTestName"].isin(aspergillus_tests)].copy()

	# Modify `mgcTestName`
	df_aspergillus["mgcTestName"] = "Aspergillus_pooled"

	# Append to original DataFrame
	df_combined = pd.concat([df, df_aspergillus], ignore_index=True)

	return df_combined

	# Define mapping function using substring rules
def map_test_name(name):
    name_lower = name.lower() if isinstance(name, str) else ''
    if 'flavus' in name_lower:
        return "Aspergillus_flavus"
    elif 'fumigatus' in name_lower:
        return "Aspergillus_fumigatus"
    elif 'niger' in name_lower:
        return "Aspergillus_niger"
    elif 'terreus' in name_lower:
        return "Aspergillus_terreus"
    elif 'salmonella' in name_lower:
        return "Salmonella"
    elif 'aerobic' in name_lower:
        return "TAC"
    elif 'yeast' in name_lower:
        return "TYM"
    elif 'aspergillus' in name_lower:
        return "Aspergillus_spp"
    elif 'stec' in name_lower or 'shiga' in name:
        return "STEC"
    elif 'Pesticides (pass/fail)' == name:
    	return "PesticidesPassFail" 
    else:
        return None

def merge_new_tests(df_existing, df_expanded):
    """
    Merge new test records into existing dataset, avoiding duplicates
    based on (SampleId, TestName).

    Parameters:
    - df_existing (pd.DataFrame): existing data with test results
    - df_expanded (pd.DataFrame): new potential tests to add (may include duplicates)

    Returns:
    - pd.DataFrame: combined DataFrame with no duplicate (SampleId, TestName)
    """

    # Ensure consistent dtypes
    df_existing = df_existing.copy()
    df_expanded = df_expanded.copy()

    # Filter rows from df_expanded that don't already exist in df_existing
    new_rows_mask = ~df_expanded.set_index(['SampleId', 'TestName']).index.isin(
        df_existing.set_index(['SampleId', 'TestName']).index
    )
    df_new_rows = df_expanded[new_rows_mask]

    # Combine and drop exact duplicates if needed
    df_merged = pd.concat([df_existing, df_new_rows], ignore_index=True)

    # Optional: drop duplicates again just to be safe (can be omitted if not needed)
    df_merged = df_merged.drop_duplicates(subset=['SampleId', 'TestName'])

    return df_merged.reset_index(drop=True)

class AutoVivification(dict):
	"""Implementation of perl's autovivification feature."""
	def __getitem__(self, item):
		try:
			return dict.__getitem__(self, item)
		except KeyError:
			value = self[item] = type(self)()
			return value

def load_file_drop_nulls(input_csv):
	
	input_csv = '/Dragen/CCC_Lab_Testing_Data/OR/Anonymized Test Data Feb 2021 to April 2024.csv'

	df = pd.read_csv(input_csv)
	df = df.dropna(how="all") #if all values are null drop (i.e. a blank line)
	df = df.replace(r'^\s*$', np.nan, regex=True)
	
	return df

def add_date_tested(df):

	df['TestPerformedDate'] = df['MonthOfTest'].astype(str) + '/1/' + df['YearOfTest'].astype(str)
	df['TestPerformedDate'] = pd.to_datetime(df['TestPerformedDate'])
	df['TestPerformedDate'].dt.strftime('%Y-%m-%d')

	return df

def add_mgc_test_name(df):
	# Apply mapping
	df['mgcTestName'] = df['TestName'].apply(map_test_name)
	return df

def remove_RD_samples(df):
	# Example regex pattern matching 'R&D' anywhere in TestName
	pattern = r'R&D'

	# Select SampleIds where TestName matches the regex
	sample_ids_to_exclude = df.loc[df['TestName'].str.contains(pattern, regex=True, na=False), 'SampleId'].unique()

	# Remove all rows with those SampleIds
	df_filtered = df[~df['SampleId'].isin(sample_ids_to_exclude)]
	return df_filtered

def make_df_all_pass_sample_test_combos(df):

	df_distinct_sampleID = (
	    df.groupby('SampleId')
	    .agg({
	        'LabId': lambda x: int(x.iloc[0]),
	        'SampleCreatedByLicenseId': lambda x: x.iloc[0],
	        'MonthOfTest': lambda x: x.iloc[0],
	        'YearOfTest': lambda x: x.iloc[0],
	        'ProductType': lambda x: x.iloc[0],
	        'TestPerformedDate': lambda x: x.iloc[0]
	    })
	    .reset_index()
	)

	distinct_test_names = df['TestName'].dropna().unique()
	test_names_df = pd.DataFrame({'TestName': distinct_test_names})

	# Add key for cross join
	df_distinct_sampleID['_key'] = 1
	test_names_df['_key'] = 1

	# Cross join: every SampleId x every TestName
	df_expanded = pd.merge(df_distinct_sampleID, test_names_df, on='_key').drop('_key', axis=1)

	# Add empty/missing Result and PassFail columns
	df_expanded['Result'] = pd.NA #when creating a test entry w/a PASS for test data we don't have, there is no Result so we set to NULL
	df_expanded['PassFail'] = True

	# Optional: reorder columns if desired
	cols_order = ['SampleId', 'LabId', 'SampleCreatedByLicenseId', 'MonthOfTest', 'YearOfTest',
	              'ProductType', 'TestPerformedDate', 'TestName', 'Result', 'PassFail']

	df_expanded = df_expanded[cols_order]
	return df_expanded

def cleanup_dataframe_for_mysql_load_table(df):
	# Clean up: drop blank rows and replace blank/whitespace strings
	df = df.dropna(how="all")  # Drop rows where all values are NaN
	df = df.replace(r'^\s*$', np.nan, regex=True)
	df = df.fillna('\\N')  # Replace NaNs with string '\N'

	# Lastly, get rid of everything that doesn't have an mgcTestName and won't get dumped into the labtesting portal SQL table:
	df = df[df['mgcTestName'] != '\\N']

	# and finally, we need to make it look like our SQL table:

	#mysql> desc states_merged;
	#+-------------------+--------------+------+-----+---------+-------+
	#| Field             | Type         | Null | Key | Default | Extra |
	#+-------------------+--------------+------+-----+---------+-------+
	#| LabName           | varchar(250) | YES  | MUL | NULL    |       |
	#| State             | varchar(2)   | YES  | MUL | NULL    |       |
	#| TestPassed        | tinyint(1)   | YES  |     | NULL    |       |
	#| TestResultLevel   | float        | YES  |     | NULL    |       |
	#| TestPerformedDate | date         | YES  |     | NULL    |       |
	#| StateTestName     | varchar(250) | YES  |     | NULL    |       |
	#| mgcTestName       | varchar(250) | YES  | MUL | NULL    |       |
	#| StateMaterialType | varchar(250) | YES  |     | NULL    |       |
	#+-------------------+--------------+------+-----+---------+-------+
	#8 rows in set (0.01 sec)

	df['State'] = 'OR'

	df = df.rename(columns={
	    'LabId': 'LabName',  # if LabId is used as LabName
	    'Result': 'TestResultLevel',
	    'PassFail': 'TestPassed',
	    'TestName': 'StateTestName',
	    'ProductType': 'StateMaterialType',
	})

	df['TestPassed'] = df['TestPassed'].replace({'True': 1, 'False': 0, '\\N': None}).astype('float').astype('Int64')
	df['TestResultLevel'] = pd.to_numeric(df['TestResultLevel'], errors='coerce')
	df['TestPerformedDate'] = pd.to_datetime(df['TestPerformedDate'], errors='coerce').dt.date

	df = df[['LabName', 'State', 'TestPassed', 'TestResultLevel',
	         'TestPerformedDate', 'StateTestName', 'mgcTestName', 'StateMaterialType']]

	final_df = create_aspergillus_pooled_category(df)

	final_df = final_df.dropna(how="all")  # Drop rows where all values are NaN
	final_df = final_df.replace(r'^\s*$', np.nan, regex=True)
	final_df = final_df.fillna('\\N')  # Replace NaNs with string '\N'

	return final_df


def main():
	#$ md5sum "/Dragen/CCC_Lab_Testing_Data/OR/Anonymized Test Data Feb 2021 to April 2024.csv"
	#3ebb15bfcee04cf35092852d00333edd  /Dragen/CCC_Lab_Testing_Data/OR/Anonymized Test Data Feb 2021 to April 2024.csv
	input_csv = '/Dragen/CCC_Lab_Testing_Data/OR/Anonymized Test Data Feb 2021 to April 2024.csv'


	print("Load File and Drop nulls...", file=sys.stderr)
	df = load_file_drop_nulls(input_csv=input_csv)
	print("Add date tested in TimeDate format...", file=sys.stderr)
	df = add_date_tested(df=df)
	print("Add MGC TestNames used for standardizing test names and comparing states...", file=sys.stderr)
	df = add_mgc_test_name(df=df)
	print("Remove any tests matching 'R&D'...", file=sys.stderr)
	df = remove_RD_samples(df=df)

	#next, the logic being implemented is:
	# All SampleIDs that: 
	# 1. exist
	# 2. aren't R&D samples (which have already been removed)
	# 3. are missing any tests
	# 4. SHOULD actually have PassFail results of "True" indicating that the passed the test that is missing

	#first, make a dataframe that has every SampleID/TestName prepopulated with a True (PassFail value) for every test:

	print("Make a dummy dataframe with a PASS for all SampleID/TestName combos...", file=sys.stderr)
	df_all_pass = make_df_all_pass_sample_test_combos(df=df)

	# If a SampleID/TestName combo does not exist in df, add a passing test with Null values for "Result" (as we have no way of knowing what that is).
	# This should be okay as the labtesting portal uses True/False in the PassFail field and ignores Result:

	print("Merge the fake PASS data into the real dataframe to create PASS test results for missing data only...", file=sys.stderr)
	df = merge_new_tests(df_existing=df, df_expanded=df_all_pass)

	print("Adding mgcTestNames for new data missing this...", file=sys.stderr)
	# now we have new values that are missing an mgcTestName so let's add those again:
	df['mgcTestName'] = df['TestName'].apply(map_test_name)

	#print("Write dataframe to tab-delimited file...\n", file=sys.stderr)
	#print(df.to_csv(sep='\t', index=False, header=True))
	#exit()

	print("Cleaning up dataframe to make it play nice with 'LOAD DATA LOCAL INFILE' so it can be added to 'states_merged' SQL table...", file=sys.stderr)
	final_df = cleanup_dataframe_for_mysql_load_table(df=df)

	print(final_df.to_csv(sep='\t', index=False, header=True))

	exit()

if __name__ == "__main__":
	main()

