# %%
import os
import requests
import pandas as pd
from sodapy import Socrata

os.chdir('c:\\Users\\John\\Documents\\GitHub')
import keyring_login

os.chdir('c:\\Users\\John\\Documents\\GitHub\\crimeBQR')

data_dict = {
    'arrests_hist': [0, '8h9b-rp9u'],
    'arrests_ytd': [1, 'uip8-fykc'],
    'complaints_hist': [2, 'qgea-i56i'],
    'complaints_ytd': [3, '5uac-w243'],
    'mvc_crashes': [4, 'h9gi-nx95'],
    'mvc_vehicles': [5, 'bm4k-52h4'],
    'mvc_person': [6, 'f55k-p6yu'],
    'shootings_hist': [7, '833y-fsy8'],
    'shootings_ytd': [8, '5ucz-vwe8'],
    'court_summons_hist': [9, 'sv2w-rv3k'],
    'court_summons_ytd': [10, 'mv4k-y93f']
}
dict_list = list(data_dict.items())

def auth_socrata(api_type = None):

    if api_type == None:
        print(dict_list)
        selection = input('Please select a number corresponding to a crime type')
        api_type = dict_list[int(selection)][1][1]

    login_cred = keyring_login.login()

    client = Socrata("data.cityofnewyork.us",
                str(login_cred[2]),
                username = login_cred[0],
                password = login_cred[1]
    )

    # Results returned as JSON from API / converted to Python list of dictionaries by sodapy.
    results = client.get(api_type, limit=500)#5012956)
    results_df = pd.DataFrame.from_records(results)

    return(results_df)
# %%

urls = {}
for crime in data_dict.values():
    urls[crime] = requests.get('https://data.cityofnewyork.us/resource/' + str(crime) + '.json')

api_type = []
for crime in data_dict.values():
    api_type[crime] = crime
