import pandas as pd
import requests
import re
import json


def find_drugbank_id(data):
    # Regular expression to match 'https://www.drugbank.ca/drugs/' followed by 'DB' and digits
    pattern = re.compile(r'https://www\.drugbank\.ca/drugs/(DB\d+)')
    
    if isinstance(data, dict):
        for key, value in data.items():
            if isinstance(value, str):
                match = pattern.search(value)
                if match:
                    print(f"match: {match}")
                    return match.group(1)  # Return just the DB ID part
            """elif isinstance(value, (dict, list)):
                
                result = find_drugbank_id(value)
                if result:
                    return result
    elif isinstance(data, list):
        print("data is a list")
        for item in data:
            result = find_drugbank_id(item)
            if result:
                return result"""
    return None

def cid_to_drugbank(cid):

    #this is the pattern of the web link. from cid it wil be linked to the pubchem page, and in the content of the page we wanna find drugbank ID
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        data = response.json()
        
        # Search for DrugBank ID in the entire JSON data
        drugbank_id = find_drugbank_id(data)
        print(drugbank_id)
        return drugbank_id
    except:
        return None



def cid_to_drugbankid(druginfo):

    #We just want to search for drugbank IDs, where the Value of which in drugbank_id column is reported null but not null for cid, cause we need it to convert it to drugbankID
    nulldrugbank = druginfo[druginfo['drugbank_id'].isna()]
    nulldrugbank = nulldrugbank[pd.notna(nulldrugbank['cid'])]

    #Want to chase the source where we fetch drugbank ID
    druginfo['drugbankID_Source_pubchem'] = None


    #want to save the number of drugs dont have drugbank ID to compare after fetching data
    initial_null = len(nulldrugbank)


    for index, row in nulldrugbank.iterrows():
        
        cid = row['cid']
        # cid could be saved as a float datatype
        try:
            cid = int(cid)
        except ValueError:
            print(f"Value {cid} cannot be converted to an integer, continuing...")
            continue

        drugbank_id = cid_to_drugbank(cid)
        
        #if drugbank_id had a value, not "None"
        if drugbank_id:

            #fill the drugbank ID with the corresponding ID
            druginfo.loc[index, "drugbank_id"] = drugbank_id
            
            # want to store drugbank ID had come from which source
            druginfo.loc[index, "drugbankID_Source_pubchem"] = True

    return druginfo, initial_null, druginfo['drugbank_id'].isna().sum()

def anomaly(data):
    
    if data.loc[5695,'dname'] == 'myo-inositol':
        data.loc[5695,'drugbank_id'] = 'DB13178'

    if data.loc[5850,'dname'] == '32462-30-9':
        data.loc[5850,'drugbank_id'] = 'DB04291'

    if data.loc[6465,'dname'] == '1y6a':
        data.loc[6465,'drugbank_id'] = 'DB07334'

    return data


def unichem_inchikey_to_dbID(druginfo):

    
    druginfo['drugbankID_Source_inchikeys(unichem)'] = None


    di_null_drugbank = druginfo[druginfo["drugbank_id"].isna()]
    di_null_db_notnullInchikeys=di_null_drugbank[pd.notna(di_null_drugbank["inchikey"])]

    for i,row in di_null_db_notnullInchikeys.iterrows():
    
        inchikey = row['inchikey']

        # Construct the URL from inchikey
        base_url = "https://www.ebi.ac.uk/unichem/rest/verbose_inchikey"
        url = f"{base_url}/{inchikey}"
        
        response = requests.get(url)
        
        try:
            # Parse the JSON response
            data = response.json()
            
            # Look for DrugBank ID
            drugbank_id = None
            for source in data:
                if 'src_url' in source:
                    match = re.search(r'http://www\.drugbank\.ca/drugs/(DB\d+)', source['src_url'])
                    if match:
                        drugbank_id = match.group(1)
                        break
            
            if drugbank_id:
                #print(f"DrugBank ID: {drugbank_id}")
                if druginfo.loc[i, "inchikey"] == inchikey:
                    druginfo.loc[i, 'drugbank_id'] = drugbank_id
                    druginfo.loc[i, 'drugbankID_Source_inchikeys(unichem)'] = True

        except:
            print("e")
                    
    return druginfo, druginfo['drugbank_id'].isna().sum()


def unichem_drugbank(compund, sourceID):
    
    url = "https://www.ebi.ac.uk/unichem/api/v1/connectivity"

    payload = {
        "compound": compund ,# CHEMBLE ID  EXAMPLE: "CHEMBL53463"
        "searchComponents": False,
        "sourceID": sourceID,
        "type": "sourceID"
    }
    
    headers = {
        "Content-Type": "application/json"
    }
    
    response = requests.post(url, json=payload, headers=headers)
    
    if response.status_code == 200:
        data = response.json()
        
        # Convert the data to a string
        data_str = json.dumps(data)
        
        # Use regex to find the DrugBank ID
        pattern = r'"url": "http://www\.drugbank\.ca/drugs/(DB\d+)"'
        match = re.search(pattern, data_str)
        
        if match:
            drugbank_id = match.group(1)
            return drugbank_id
        else:
            return None
    else:
        return None

    

def unichem_chemblID_to_dbID(druginfo):

    druginfo['drugbankID_Source_chembl(unichem)'] = None
    
    #only fill those rows having null value in drugbank_id
    null_drugbank = druginfo[(druginfo['drugbank_id'].isna())]

    # As it wants to convert chembl to drugbank, you must work with the subset of the data not having null value in chemblID
    null_drugbank = null_drugbank[pd.notna(null_drugbank['chembl_id'])]

    for i, row in null_drugbank.iterrows():
        
        chembl = row['chembl_id']

        #source ID for chembl in unichem is 1
        drugbankid = unichem_drugbank(chembl, "1")
        
        if drugbankid is not None:
            
            if druginfo.loc[i, "chembl_id"] == chembl:
                druginfo.loc[i, "drugbank_id"] = drugbankid

                # To chase the source of the data in that cell
                druginfo.loc[i, "drugbankID_Source_chembl(unichem)"] = True

    return druginfo, druginfo['drugbank_id'].isna().sum()


def unichem_keggID_to_dbID(druginfo):
    
    druginfo['drugbankID_Source_kegg(unichem)'] = None
    
    #only fill those rows having null value in drugbank_id
    null_drugbank = druginfo[(druginfo['drugbank_id'].isna())]

    # As it wants to convert chembl to drugbank, you must work with the subset of the data not having null value in chemblID
    null_drugbank = null_drugbank[pd.notna(null_drugbank['kegg_id'])]

    for i, row in null_drugbank.iterrows():
        kegg = row['kegg_id']
        
        if pd.notna(kegg):

            #source ID for chembl in unichem is 6
            drugbankid = unichem_drugbank(kegg, "6")
            
            if drugbankid is not None:
                
                if druginfo.loc[i, "kegg_id"] == kegg:
                    druginfo.loc[i, "drugbank_id"] = drugbankid
                    druginfo.loc[i, "drugbankID_Source_kegg(unichem)"] = True

    return druginfo, druginfo['drugbank_id'].isna().sum()



# Load your DataFrame- This is the druginfo table from drugcomb api for all the existing drugs having drugnames, inchikeys, 
druginfo = pd.read_csv("/g/data/yr31/pm5363/2024/Codes/Drugs_api.csv")



druginfo, initial_null, null_dbID_after_convert_cid_to_dbID = cid_to_drugbankid(druginfo)


#Some samples had drugbankIDs that dont match the pattern like: 'DB03106; DB13178' output for only one sample, I manually searched wich one is the correct dbID for that drug and updated it
druginfo = anomaly(druginfo)

druginfo, null_dbID_after_convert_inchikeys_to_dbID = unichem_inchikey_to_dbID(druginfo)

                              
druginfo, null_dbID_after_convert_chembl_to_dbID = unichem_chemblID_to_dbID(druginfo)


druginfo, null_dbID_after_convert_kegg_to_dbID = unichem_keggID_to_dbID(druginfo)

print(f"initial null: {initial_null}\n null_dbID_after_convert_cid_to_dbID: {null_dbID_after_convert_cid_to_dbID}\n null_dbID_after_convert_inchikeys_to_dbID: {null_dbID_after_convert_inchikeys_to_dbID}\n null_dbID_after_convert_chembl_to_dbID: {null_dbID_after_convert_chembl_to_dbID}\n null_dbID_after_convert_kegg_to_dbID: {null_dbID_after_convert_kegg_to_dbID}")
druginfo.to_csv("/g/data/yr31/pm5363/2024/Codes/Drugs_fill_na.csv", index=False)
                              

