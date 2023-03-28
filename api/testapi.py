import requests

url = 'http://localhost:5000/standardize'

response = requests.post(url, json={'smiles_list': []})
print (response.json())

long_smiles = "".join(['C' for i in range(10)])
response = requests.post(url, json={'smiles_list': [{"smiles": long_smiles}]})
print(response.json())

invalid_list = [
    {'smiles':'invalid','id':'1'},
    {'smiles':'SMILES','id':'2'},
    {'smiles':'strings','id':'3'}
]
response = requests.post(url, json={'smiles_list': invalid_list})
print(response.json())  # Should return an error message or an empty result

duplicate_list = [
    {"smiles": "C", "id": "1"},
    {"smiles": "CC", "id": "2"},
    {"smiles": "C", "id": "3"},
    {"smiles": "CCC", "id": "4"}
]
response = requests.post(url, json={'smiles_list': duplicate_list})
print(response.json())

aminophylline = [{"smiles":"Cn1c2nc[nH]c2c(=O)n(C)c1=O", "id": "aminophylline"}]
response = requests.post(url, json={'smiles_list': aminophylline})
print(response.json())
